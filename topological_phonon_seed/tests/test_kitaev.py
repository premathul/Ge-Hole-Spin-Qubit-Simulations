import numpy as np
from topological_phonon_seed.kitaev_strain import (bdg_chain, particle_hole_residual,
                                                   majorana_number, bulk_gap, control_cycle, diagnostics)

def test_phs_hermiticity_and_spectrum():
    for theta in (0.,.7,2.):
        H=bdg_chain(12,.3,pairing=.6*np.exp(1j*theta))
        assert np.linalg.norm(H-H.conj().T)<1e-12
        assert particle_hole_residual(H)<1e-12
        e=np.linalg.eigvalsh(H)
        assert np.allclose(e,-e[::-1],atol=1e-11)

def test_static_invariant_and_gap_closing():
    assert majorana_number(0)==-1
    assert majorana_number(3)==1
    assert majorana_number(2)==0
    assert bulk_gap(0)>0
    assert bulk_gap(2)<1e-10

def test_no_false_pumping_claim():
    mu,d=control_cycle(np.linspace(0,2*np.pi,101))
    assert np.isclose(mu[0],mu[-1]) and np.isclose(d[0],d[-1])
    result=diagnostics(mu,d)
    assert result["minimum_bulk_gap"]>0
    assert result["parity_pump_invariant"]=="NOT COMPUTED"
