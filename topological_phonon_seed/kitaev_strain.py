"""Toy spinless Kitaev BdG model. NOT microscopic Ge or a parity-pump calculation."""
import numpy as np

def bdg_chain(n, mu, hopping=1., pairing=0.6):
    if n < 2: raise ValueError("n >= 2 required")
    h = -mu*np.eye(n, dtype=complex)
    d = np.zeros((n,n), dtype=complex)
    for j in range(n-1):
        h[j,j+1] = h[j+1,j] = -hopping
        d[j,j+1], d[j+1,j] = pairing, -pairing
    return np.block([[h,d],[d.conj().T,-h.T]])

def particle_hole_residual(H):
    n = H.shape[0]//2
    x = np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])
    return float(np.linalg.norm(x @ H.conj() @ x + H))

def majorana_number(mu, hopping=1.):
    """Static class-D number: -1 topological, +1 trivial, 0 at transition."""
    return int(np.sign((-mu-2*hopping)*(-mu+2*hopping)))

def bulk_gap(mu, hopping=1., pairing=0.6, nk=2048):
    k = np.linspace(-np.pi,np.pi,nk,endpoint=False)
    return float(np.sqrt((-mu-2*hopping*np.cos(k))**2
                          +4*abs(pairing)**2*np.sin(k)**2).min())

def control_cycle(t, mu0=.3, delta0=.6, a_mu=.2, a_delta=.1, phase=np.pi/2):
    """Phenomenological controls; NOT derived from Ge Bir-Pikus coupling."""
    return mu0+a_mu*np.cos(t), delta0+a_delta*np.cos(t+phase)

def diagnostics(mu, delta):
    gaps = [bulk_gap(float(m), pairing=float(d)) for m,d in zip(mu,delta)]
    return {"minimum_bulk_gap":float(min(gaps)),
            "static_invariants":sorted(set(majorana_number(float(m)) for m in mu)),
            "parity_pump_invariant":"NOT COMPUTED"}

if __name__ == "__main__":
    import json
    m,d = control_cycle(np.linspace(0,2*np.pi,101))
    print(json.dumps(diagnostics(m,d),indent=2))
