"""Resumable, overlap-tracked low-field Zeeman study on existing meshes.

Recovered from the original AnnularArc_SweetSpot_DQD_Ge project record.
"""
from pathlib import Path
import gc, json, os, sys, time
import numpy as np
from qtcad.device import Device, SubDevice, constants as ct, materials as mt, io, analysis as an
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device.schrodinger import Solver as SchrodingerSolver, SolverParams

W = Path(__file__).resolve().parents[1]
O = W / "33_zeeman_linear_tracking"
O.mkdir(exist_ok=True)
sys.path.insert(0, str(W / "parameters"))
import params as p

GATES = ["RES_L", "LB", "L", "CB", "R", "RB", "RES_R"]
DOTS = ["SiGe_bottom.dot_region","Ge_well.dot_region","SiGe_top.dot_region","Ge_cap.dot_region"]
DIRS = {
    "x":[1,0,0],"y":[0,1,0],"z":[0,0,1],
    "xy":[1,1,0],"xz":[1,0,1],"yz":[0,1,1],
    "xyz":[1,1,1],"xmy":[1,-1,0],"123":[1,2,3],
}
DIRS = {k:(np.asarray(v,float)/np.linalg.norm(v)).tolist() for k,v in DIRS.items()}
FIELDS_MT = [0.5,1.0,2.0,5.0,10.0,20.0]
MESHES = ["coarse","baseline","fine"]
NSTATES = 12
TOL_EV = 1e-12

def build(mesh_name):
    d = Device(
        Mesh(1e-9, W/"meshes"/f"annulararc_dqd_ge_{mesh_name}.msh"),
        conf_carriers="h",
        hole_kp_model="luttinger_kohn_foreman",
    )
    d.set_temperature(p.TEMPERATURE_K)
    sige = mt.Material(alloy_func=mt.gen_SiGe_DFT, alloy_composition=p.X_GE_BARRIER)
    sige.hole_kp_params = mt.alloy_bowing("hole_kp_params", mt.Ge, mt.Si, p.X_GE_BARRIER)
    sige.hole_Zeeman_params = mt.alloy_bowing("hole_Zeeman_params", mt.Ge, mt.Si, p.X_GE_BARRIER)
    sige.vlnce_band_def_pot = mt.alloy_bowing("vlnce_band_def_pot", mt.Ge, mt.Si, p.X_GE_BARRIER)
    ge = mt.Material(inp_dict=mt.gen_Ge_strained_on_SiGe(p.X_GE_BARRIER))
    ge.hole_kp_params = p.HOLE_KP_PARAMS
    ge.hole_Zeeman_params = np.array([p.HOLE_ZEEMAN_KAPPA,p.HOLE_ZEEMAN_Q])
    ge.vlnce_band_def_pot = np.array([p.DEFORM_A_V_EV,p.DEFORM_B_V_EV,p.DEFORM_D_V_EV])*ct.e
    for name,mat in [
        ("substrate",sige),("SiGe_bottom",sige),("Ge_well",ge),("SiGe_top",sige),
        ("Ge_cap",mt.Ge),("high-k_gate",mt.Al2O3),
        ("SiGe_bottom.dot_region",sige),("Ge_well.dot_region",ge),
        ("SiGe_top.dot_region",sige),("Ge_cap.dot_region",mt.Ge)
    ]:
        d.new_region(name,mat)
    d.align_bands(ge)
    d.set_strain(p.STRAIN_TENSOR, region="Ge_well")
    d.set_strain(p.STRAIN_TENSOR, region="Ge_well.dot_region")
    for g in GATES:
        d.new_gate_bnd(g,p.GATE_VOLTAGES[g],p.GATE_WORK_FUNCTION_EV*ct.e)
    d.set_dot_region(DOTS)
    d.phi = io.load(W/"26_mesh_convergence"/mesh_name/"electrostatic_potential.hdf5", var_name="var")
    d.set_V_from_phi()
    return d

def overlap(mesh,a,b):
    q = np.sum(np.conj(a)*b,axis=1)
    return an.integrate(mesh,np.real(q)) + 1j*an.integrate(mesh,np.imag(q))

def solve(mesh_name,B):
    d = build(mesh_name)
    d.set_Bfield(np.asarray(B,float))
    d.set_Zeeman(True)
    d.set_orb_B_effects(True)
    sub = SubDevice(d,SubMesh(d.mesh,DOTS))
    sp = SolverParams()
    sp.num_states = NSTATES
    sp.tol = TOL_EV
    sp.maxiter = 2000
    sp.method = "fast"
    sp.guess = "box"
    t = time.time()
    SchrodingerSolver(sub,solver_params=sp).solve()
    elapsed = time.time()-t
    psi = np.asarray(sub.eigenfunctions)
    E = np.asarray(sub.energies)
    xyz = np.asarray(sub.mesh.xyz()).T
    left = p.node_arc_length_nm(xyz[:,0],xyz[:,1]) < 0
    dens = np.sum(abs(psi)**2,axis=2)
    norm = np.array([an.integrate(sub.mesh,dens[:,i]).real for i in range(NSTATES)])
    PL = np.array([an.integrate(sub.mesh,dens[:,i]*left).real/norm[i] for i in range(NSTATES)])
    return sub,E,psi,PL,elapsed

for mesh_name in MESHES:
    md = O/mesh_name
    md.mkdir(exist_ok=True)
    refp = md/"B0_reference.npz"
    if refp.exists():
        with np.load(refp) as q:
            E0=q["E_J"]; ref=q["ref_psi"]; ref_idx=q["ref_idx"]; PL0=q["PL"]
    else:
        sub,E0,psi0,PL0,elapsed = solve(mesh_name,[0,0,0])
        cand = np.where(PL0>.90)[0]
        ref_idx = cand[np.argsort(E0[cand])[:2]]
        if len(ref_idx)!=2:
            raise RuntimeError(f"{mesh_name}: cannot identify B=0 left-dot Kramers pair; PL={PL0}")
        ref = psi0[:,ref_idx,:].copy()
        np.savez_compressed(refp,E_J=E0,ref_psi=ref,ref_idx=ref_idx,PL=PL0,tol_eV=TOL_EV,num_states=NSTATES,elapsed_s=elapsed)
        del sub,psi0
        gc.collect()

    ck = md/"checkpoint.json"
    data = json.loads(ck.read_text()) if ck.exists() else {
        "mesh":mesh_name,"tol_eV":TOL_EV,"num_states":NSTATES,"points":{}
    }

    for direction,u0 in DIRS.items():
        u=np.asarray(u0)
        for BmT in FIELDS_MT:
            key=f"{direction}_{BmT:g}mT"
            if key in data["points"]:
                continue
            (O/"status.json").write_text(json.dumps({
                "state":"RUNNING","pid":os.getpid(),"mesh":mesh_name,"point":key
            },indent=2))
            sub,E,psi,PL,elapsed = solve(mesh_name,u*BmT*1e-3)
            U=np.array([[overlap(sub.mesh,ref[:,a,:],psi[:,k,:]) for k in range(NSTATES)] for a in range(2)])
            best=None
            for i in range(NSTATES):
                for j in range(i+1,NSTATES):
                    s=np.linalg.svd(U[:,[i,j]],compute_uv=False)
                    score=float(np.sum(s*s))
                    if best is None or score>best[0]:
                        best=(score,i,j,s)
            score,i,j,s=best
            pair=sorted((i,j),key=lambda k:E[k])
            split=float(abs(E[pair[1]]-E[pair[0]]))
            other=[k for k in range(NSTATES) if k not in pair]
            gap=float(min(abs(E[k]-E[j]) for k in other for j in pair))
            rec={
                "direction":direction,
                "u":u.tolist(),
                "B_T":BmT*1e-3,
                "pair":pair,
                "split_ueV":split/ct.e*1e6,
                "g_eff":split/(ct.muB*BmT*1e-3),
                "subspace_score":score,
                "singular_values":s.tolist(),
                "neighbor_weights":np.sum(abs(U)**2,axis=0).tolist(),
                "nearest_gap_ueV":gap/ct.e*1e6,
                "PL_pair":PL[pair].tolist(),
                "elapsed_s":elapsed,
            }
            np.savez_compressed(
                md/f"{key}.npz",
                E_J=E,PL=PL,overlap=U,
                **{k:np.asarray(v) for k,v in rec.items() if k not in ("direction",)}
            )
            data["points"][key]=rec
            ck.write_text(json.dumps(data,indent=2))
            del sub,psi
            gc.collect()

(O/"status.json").write_text(json.dumps({"state":"COMPLETE","pid":None},indent=2))
print("ZEEMAN_LINEAR_TRACKING_COMPLETE",flush=True)
