"""1325 author: one shared finite quantum reference and two opposite ports.

Finite-sum, embedding and rectangular-norm methods extend the frozen 1324
author construction; no old or independent scientific implementation is imported.
Default and --check are read-only; --write-new is exclusive.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent.parent
COST = (0, 1, 1, 2)  # D,L0,L1,U


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def same_result(a, b, path="result"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for k in a:
            same_result(a[k], b[k], path + "." + k)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for k, (x, y) in enumerate(zip(a, b)):
            same_result(x, y, f"{path}[{k}]")
    elif isinstance(a, float):
        assert isinstance(b, (float, int)) and not isinstance(b, bool), path
        assert np.isfinite(b) and np.isclose(a, b, rtol=2e-10, atol=2e-11), path
    else:
        assert type(a) is type(b) and a == b, path


def opnorm(A):
    return float(np.sqrt(max(0., np.linalg.eigvalsh(A.conj().T @ A)[-1])))


def unitary(H, time):
    w, q = np.linalg.eigh(H)
    return (q * np.exp(-1j * time * w)) @ q.conj().T


def moments(L, nu, phase=0.):
    k = np.arange(L + 1, dtype=float)
    logp = k * np.log(nu) - np.array([math.lgamma(j + 1) for j in range(L + 1)])
    logp -= logp.max()
    p = np.exp(logp)
    p /= p.sum()
    s = float(np.sqrt(nu) * np.sum(p[:-1] / np.sqrt(k[:-1] + 1)))
    z = float(np.sqrt(nu) * np.sum(p[:-1] / np.sqrt(k[:-1] + 2)))
    c = float(np.sum(p / (k + 1)))
    kap = float(nu * np.sum(p[:-2] / (k[:-2] + 2))) if L >= 2 else 0.
    return dict(p=p, coeff=np.sqrt(p) * np.exp(1j * phase * k),
                s=s, z=z, c=c, kappa=kap, mean=float(p @ k),
                top=float(p[-1]), vacuum=float(p[0]))


def occupations(Kmax, n, Kmin=0):
    basis = []
    for K in range(Kmin, Kmax + 1):
        for m in range(4):
            N = K - COST[m]
            if N >= 0:
                for n0 in range(N + 1):
                    for f in range(n):
                        basis.append((m, n0, N - n0, f))
    return basis


def parent_matrix(basis, MD, MU):
    """Direct one-fermion compression of the even CAR vertices in proof (4).

    A right annihilator meets the only occupied fermion; the left creation
    then gives coefficient +1. The epsilon sign below is independent of CAR.
    """
    pos = {x: j for j, x in enumerate(basis)}
    n = len(MD)
    V = np.zeros((len(basis), len(basis)), complex)
    for (m, n0, n1, f), j in pos.items():
        N = n0 + n1
        if m == 0 and N:
            for alpha, occ in enumerate((n0, n1)):
                if occ:
                    for h in range(n):
                        i = pos[(alpha + 1, n0 - int(alpha == 0),
                                 n1 - int(alpha == 1), h)]
                        v = np.sqrt(occ / N) * MD[h, f]
                        V[i, j] += v
                        V[j, i] += v.conjugate()
        elif m == 3:
            for alpha, beta, sign in ((0, 1, 1), (1, 0, -1)):
                occ = (n0, n1)[beta]
                for h in range(n):
                    i = pos[(alpha + 1, n0 + int(beta == 0),
                             n1 + int(beta == 1), h)]
                    v = sign * np.sqrt((occ + 1) / (N + 2)) * MU[h, f]
                    V[i, j] += v
                    V[j, i] += v.conjugate()
    return V


def source_diagonals(basis):
    nb = np.array([x[1] + x[2] for x in basis], float)
    nd = np.array([float(x[0] == 0) for x in basis])
    nl = np.array([float(x[0] in (1, 2)) for x in basis])
    nu = np.array([float(x[0] == 3) for x in basis])
    K = nb + nl + 2 * nu
    return dict(nb=nb, nd=nd, nl=nl, nu=nu, K=K,
                h0=nd + 2 * nl + 3 * nu + nb,
                charge=.23 * nd + .60 * nl + .97 * nu + .37 * nb)


def embedding(basis, n, data, u):
    J = np.zeros((len(basis), 4 * n), complex)
    L = len(data["p"]) - 1
    for row, (m, n0, n1, f) in enumerate(basis):
        N = n0 + n1
        if N <= L:
            J[row, m * n + f] = (data["coeff"][N] * np.sqrt(math.comb(N, n0))
                                 * u[0] ** n0 * u[1] ** n1)
    return J


def effective(MD, MU, data, u, phase):
    n = len(MD)
    h = np.zeros((4 * n, 4 * n), complex)
    dual = np.array([u[1].conjugate(), -u[0].conjugate()])
    for a in range(2):
        left = slice((a + 1) * n, (a + 2) * n)
        h[left, :n] = np.exp(1j * phase) * data["s"] * u[a] * MD
        h[left, 3*n:] = np.exp(-1j * phase) * data["z"] * dual[a] * MU
        h[:n, left] = h[left, :n].conj().T
        h[3*n:, left] = h[left, 3*n:].conj().T
    return h


def residual_gram(MD, MU, data, u):
    n = len(MD)
    s, z, c = data["s"], data["z"], data["c"]
    AD, AU = MD @ MD.conj().T, MU @ MU.conj().T
    G = np.zeros((4*n, 4*n), complex)
    G[:n, :n] = (1 - data["vacuum"] - s*s) * MD.conj().T @ MD
    G[3*n:, 3*n:] = (1-z*z) * MU.conj().T @ MU
    P = np.outer(u, u.conj())
    G[n:3*n, n:3*n] = (np.kron(P, (1-s*s)*AD)
                         + np.kron(np.eye(2)-P, c*AD+(1-c-z*z)*AU))
    return G


def direction_maps(K):
    W = np.zeros((2*K, K+1))
    for j in range(K+1):
        if j:
            W[j-1, j] = np.sqrt(j/K)
        if j < K:
            W[K+j, j] = np.sqrt((K-j)/K)
    J = np.zeros((2*K, max(0, K-1)))
    for j in range(K-1):
        J[j, j] = np.sqrt((K-1-j)/K)
        J[K+j+1, j] = -np.sqrt((j+1)/K)
    return W, J


class AxisParent:
    """Exact invariant hull of every matter/flavor input with u=e0.

    Per K: (D(K,0),L0(K-1,0)) and
    (D(K-1,1),L0(K-2,1),L1(K-1,0),U(K-2,0)).
    K=0/1 and the top K=L+2 are retained. All blocks derive from the vertices.
    """
    def __init__(self, L, MD, MU):
        self.L, self.n = L, len(MD)
        self.basis, self.blocks = [], []
        raw = [[(0, 0, 0)]]
        for K in range(1, L+3):
            raw.append([(0,K,0), (1,K-1,0)])
            raw.append([(0,K-1,1), (2,K-1,0)] if K == 1 else
                       [(0,K-1,1), (1,K-2,1), (2,K-1,0), (3,K-2,0)])
        for short in raw:
            rows = [(m,a,b,f) for m,a,b in short for f in range(self.n)]
            start = len(self.basis)
            self.basis.extend(rows)
            V = parent_matrix(rows, MD, MU)
            vals, vecs = np.linalg.eigh(V)
            self.blocks.append((slice(start,len(self.basis)), V, vals, vecs))
        self.pos = {x:i for i,x in enumerate(self.basis)}
        assert len(self.pos) == len(self.basis)
        self.sources = source_diagonals(self.basis)

    def apply(self, X):
        Y = np.zeros_like(X, dtype=complex)
        for sl,V,w,q in self.blocks:
            Y[sl] = V @ X[sl]
        return Y

    def evolve(self, X, time):
        Y = np.zeros_like(X, dtype=complex)
        for sl,V,w,q in self.blocks:
            Y[sl] = (q * np.exp(-1j*time*w)) @ (q.conj().T @ X[sl])
        return Y

    def embed(self, data):
        return embedding(self.basis, self.n, data, np.array([1.,0.]))

    def coherent_readout(self, state, phase):
        cross = 0j
        for (m,a,b,f),i in self.pos.items():
            if m == 1 and (2,a,b,f) in self.pos:
                cross += state[i] * state[self.pos[(2,a,b,f)]].conjugate()
        return cross, float(-2*np.real(np.exp(-2j*phase)*cross))


def diagonal_mean(X, d):
    return X.conj().T @ (d[:,None] * X)


def compute():
    checks = []
    def close(name, a, b=0., tol=4e-9):
        error = float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        assert error <= tol, (name,error,tol)
        checks.append(dict(name=name, kind="equality", residual=error,
                           tolerance=tol, passed=True))
    def bound(name, value, upper, tol=3e-10):
        assert value <= upper+tol, (name,value,upper)
        checks.append(dict(name=name, kind="upper_bound", value=float(value),
                           upper=float(upper), tolerance=tol, passed=True))

    # Explicit direction maps and finite representation witnesses; the proof
    # supplies the all-representation exclusion, not this finite list.
    for K in (1,2,3,4):
        W,J = direction_maps(K)
        close(f"direction.K{K}.W_isometry", W.T@W, np.eye(K+1))
        close(f"direction.K{K}.J_isometry", J.T@J, np.eye(max(0,K-1)))
        close(f"direction.K{K}.cross", W.T@J)
        close(f"direction.K{K}.full_left", W@W.T+J@J.T, np.eye(2*K))
    def product_spins(a,b):
        return list(range(abs(a-b),a+b+1,2))
    role_examples=[]
    for port, reference in (("D",4),("U",2)):
        for two_j in range(7):
            spins=product_spins(reference,two_j)
            allowed=all(j in (2,4) for j in spins)
            close(f"role.K4.{port}.two_j{two_j}", int(allowed), int(two_j==0))
            role_examples.append(dict(port=port,two_j=two_j,product_spins=spins,
                                      compatible=allowed))
    low_input=product_spins(3,0)+product_spins(1,0)+product_spins(1,2)
    close("role.K3.full_coverage_counterexample", sorted(low_input), [1,1,3,3])

    rng=np.random.default_rng(1325)
    n=3
    MD=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.035
    MU=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.027
    u=np.array([np.sqrt(.39),np.sqrt(.61)*np.exp(.43j)])
    small=[]
    for L,nu in ((1,.7),(3,1.7)):
        phase=.31
        data=moments(L,nu,phase)
        basis=occupations(L+2,n)
        V=parent_matrix(basis,MD,MU)
        src=source_diagonals(basis)
        J=embedding(basis,n,data,u)
        h=effective(MD,MU,data,u,phase)
        G=residual_gram(MD,MU,data,u)
        delta=V@J-J@h
        eta=float(np.sqrt(max(0.,np.linalg.eigvalsh(G)[-1])))
        prefix=f"full.L{L}"
        close(prefix+".input_isometry",J.conj().T@J,np.eye(4*n))
        close(prefix+".Hermitian",V,V.conj().T)
        close(prefix+".mean_generator",J.conj().T@V@J,h)
        close(prefix+".general_complex_rotated_Gram",delta.conj().T@delta,G)
        close(prefix+".full_parent_norm",opnorm(V),max(opnorm(MD),opnorm(MU)))
        for key in ("K","charge","h0"):
            close(prefix+".commuting_"+key,src[key][:,None]*V,V*src[key][None,:])
        time=1.43
        X=unitary(V,time)@J
        Y=J@unitary(h,time)
        bound(prefix+".full_input_Duhamel",opnorm(X-Y),time*eta)
        close(prefix+".actual_isometry",X.conj().T@X,np.eye(4*n))
        for key in ("K","charge","h0"):
            close(prefix+".ledger_"+key,diagonal_mean(X,src[key]),diagonal_mean(J,src[key]))
        close(prefix+".interaction_energy",X.conj().T@V@X,J.conj().T@V@J)
        H=np.diag(src["h0"])+V
        close(prefix+".total_energy",X.conj().T@H@X,J.conj().T@H@J)
        close(prefix+".physical_frame",unitary(H,time)@J,
              np.exp(-1j*time*src["h0"])[:,None]*X,1e-8)
        partner=rng.normal(size=(4*n,2))+1j*rng.normal(size=(4*n,2))
        partner/=np.linalg.norm(partner)
        bound(prefix+".all_ports_unknown_partner",np.linalg.norm((X-Y)@partner),time*eta)
        # Four independent source parameters; no resonance retuning.
        for key in ("nd","nl","nu","nb"):
            eps=1e-5
            plus=np.diag(src["h0"]+eps*src[key])+V
            minus=np.diag(src["h0"]-eps*src[key])+V
            close(prefix+".independent_source_"+key,(plus-minus)/(2*eps),np.diag(src[key]),1e-8)
        # Actual reference projection and Taylor curvature, not a reversed
        # Duhamel claim about all finite-time instruments.
        tau=.037
        Z=unitary(V,tau)@J
        leak=Z-J@(J.conj().T@Z)
        gmax=max(opnorm(MD),opnorm(MU))
        bound(prefix+".finite_projection_curvature",
              abs(opnorm(leak)/tau-eta),tau*gmax*gmax/2)
        small.append(dict(L=L,nu=nu,dimension=len(basis),eta=eta,
                          full_joint_norm=opnorm(X-Y),time=time))

    # K=0/1 physical boundaries, not artificial negative boson states.
    for K in (0,1):
        basis=occupations(K,3,K)
        V=parent_matrix(basis,MD,MU)
        close(f"boundary.K{K}.no_U",sum(x[0]==3 for x in basis),0)
        if K==0:
            close("boundary.K0.complete_D_dark",V)
        else:
            close("boundary.K1.dimension",len(basis),12)

    # Full unknown-input deterministic transfer, including phases, in the
    # common odd-singular-value examples and scalar n=4/5 continuations.
    Md=np.diag([.13,.03,.01]).astype(complex)
    Mu=np.diag([.11,.05,.01]).astype(complex)
    reverse=np.eye(3)[::-1]
    transfer=[]
    common_time=np.pi/.02
    examples=[("aligned",Md,Mu),("reverse",Md,reverse@Mu),
              ("n4",.01*np.eye(4),.01*np.eye(4)),
              ("n5",.01*np.eye(5),.01*np.eye(5))]
    for label,D,U in examples:
        n=len(D);K=4
        basis=occupations(K,n,K)
        V=parent_matrix(basis,D,U)
        rights=[i for i,x in enumerate(basis) if x[0] in (0,3)]
        lefts=[i for i,x in enumerate(basis) if x[0] in (1,2)]
        X=unitary(V,common_time)[:,rights]
        T=V[np.ix_(lefts,rights)]
        W,Jdir=direction_maps(K)
        close("transfer."+label+".same_parent_direction_map",
              T,np.hstack([np.kron(W,D),np.kron(Jdir,U)]))
        a,s,vh=np.linalg.svd(T,full_matrices=False)
        expected=-1j*(a*np.sin(common_time*s))@vh
        close("transfer."+label+".no_right_failure",X[rights])
        close("transfer."+label+".all_unknown_inputs",X[lefts].conj().T@X[lefts],np.eye(len(rights)))
        close("transfer."+label+".odd_spectrum_phases",X[lefts],expected,1e-8)
        transfer.append(dict(label=label,n=n,K=K,time=float(common_time),
                             input_dimension=len(rights),failure_norm=opnorm(X[rights])))

    # Verify the axis reduction directly against the full small occupation
    # Hamiltonian before using it for the declared larger windows.
    axis=AxisParent(3,MD,MU)
    full=occupations(5,3);Vfull=parent_matrix(full,MD,MU)
    position={x:i for i,x in enumerate(full)}
    Q=np.zeros((len(full),len(axis.basis)),complex)
    for i,x in enumerate(axis.basis):Q[position[x],i]=1
    data=moments(3,1.7,.21);J=axis.embed(data)
    close("axis.full_embedding",Q@J,embedding(full,3,data,np.array([1.,0.])))
    close("axis.full_invariant_hull",Vfull@Q,Q@axis.apply(np.eye(len(axis.basis))))
    close("axis.full_evolution",Q@axis.evolve(J,1.3),unitary(Vfull,1.3)@Q@J)

    # Actual coherent cross-port readout, all final reference branches kept.
    L,nu,phase,gamma=160,64.,.27,.1
    data=moments(L,nu,phase)
    model=AxisParent(L,np.array([[gamma]]),np.array([[gamma]]))
    J=model.embed(data);time=np.pi/(2*gamma)
    X=model.evolve(J,time)
    initial=np.array([1.,0.,0.,1.])/np.sqrt(2)
    state=X@initial
    cross,readout=model.coherent_readout(state,phase)
    close("coherent.kappa_sum_identity",data["kappa"],1-data["c"]-L/(L+1)*data["top"])
    close("coherent.actual_cross_element",cross,-.5*np.exp(2j*phase)*data["kappa"])
    close("coherent.actual_readout",readout,data["kappa"])
    close("coherent.all_input_isometry",X.conj().T@X,np.eye(4))
    for key in ("K","charge","h0"):
        close("coherent.ledger_"+key,
              diagonal_mean(X,model.sources[key]),diagonal_mean(J,model.sources[key]),2e-8)
    close("coherent.interaction_energy",X.conj().T@model.apply(X),J.conj().T@model.apply(J))
    close("coherent.D_occupancy",float(diagonal_mean(X,model.sources["nl"])[0,0].real),1-data["vacuum"])
    close("coherent.U_occupancy",float(diagonal_mean(X,model.sources["nl"])[3,3].real),1.)
    gram=residual_gram(np.array([[gamma]]),np.array([[gamma]]),data,np.array([1.,0.]))
    h=effective(np.array([[gamma]]),np.array([[gamma]]),data,np.array([1.,0.]),phase)
    eta=float(np.sqrt(np.linalg.eigvalsh(gram)[-1]))
    delta=model.apply(J)-J@h
    close("coherent.all_input_Gram",delta.conj().T@delta,gram)
    bound("coherent.full_poststate_Duhamel",opnorm(X-J@unitary(h,time)),time*eta)
    # The 3-point phase twirl is exactly dephasing for this particular readout,
    # since its only contributing number difference is two. Fixed-N tests
    # additionally verify this without assuming a mean-field identity.
    phase_readouts=[]
    for shift in (0.,2*np.pi/3,4*np.pi/3):
        d=moments(L,nu,phase+shift)
        st=model.evolve(model.embed(d),time)@initial
        phase_readouts.append(model.coherent_readout(st,phase)[1])
    close("phase_average.actual_readout_zero",np.mean(phase_readouts))
    fixed=[]
    for N in (0,1,2,64):
        coeff=np.zeros(L+1,complex);coeff[N]=1
        d=dict(p=np.abs(coeff)**2,coeff=coeff)
        st=model.evolve(model.embed(d),time)@initial
        value=model.coherent_readout(st,phase)[1]
        close(f"fixed_number.N{N}.actual_readout_zero",value)
        fixed.append(dict(N=N,readout=value))
    coherent=dict(L=L,nu=nu,phase=phase,gamma=gamma,time=float(time),
                  kappa_expected=data["kappa"],actual_kappa=readout,
                  axis_dimension=len(model.basis),
                  cross_element=[float(cross.real),float(cross.imag)],
                  phase_twirl_readouts=phase_readouts,fixed_number=fixed)

    # Same reference, spectra, clock, readout and preparation; only the
    # relative left flavor orientation changes. No parameter is re-fitted.
    budget=[]
    L,nu,phase,tau,beta=64,16.,.19,.1,.0337
    data=moments(L,nu,phase)
    for label,U in (("aligned",Mu),("reverse",reverse@Mu)):
        model=AxisParent(L,Md,U)
        J=model.embed(data)
        h=effective(Md,U,data,np.array([1.,0.]),phase)
        G=residual_gram(Md,U,data,np.array([1.,0.]))
        eta=float(np.sqrt(np.linalg.eigvalsh(G)[-1]))
        delta=model.apply(J)-J@h
        close("budget."+label+".complete_Gram",delta.conj().T@delta,G)
        close("budget."+label+".same_spectra",np.linalg.svd(U,compute_uv=False),np.array([.11,.05,.01]))
        col=2*3  # L1, flavor zero is the declared common worst input.
        close("budget."+label+".specified_input_worst",float(G[col,col].real),eta*eta)
        Z=model.evolve(J,tau)
        leak=Z-J@(J.conj().T@Z)
        probability=float(np.vdot(leak[:,col],leak[:,col]).real)
        measured_rate=float(np.sqrt(probability)/tau)
        worst_rate=opnorm(leak)/tau
        remainder=tau*.13**2/2
        bound("budget."+label+".finite_time_specified_bound",abs(measured_rate-eta),remainder)
        bound("budget."+label+".finite_time_worst_bound",abs(worst_rate-eta),remainder)
        if label=="aligned":
            bound("budget.aligned.actual_rejection",beta,measured_rate)
            bound("budget.aligned.interval_rejection",beta,eta-remainder)
        else:
            bound("budget.reverse.actual_acceptance",measured_rate,beta)
            bound("budget.reverse.interval_acceptance",eta+remainder,beta)
        for key in ("K","charge","h0"):
            close("budget."+label+".ledger_"+key,
                  diagonal_mean(Z,model.sources[key]),diagonal_mean(J,model.sources[key]),2e-8)
        close("budget."+label+".interaction_energy",Z.conj().T@model.apply(Z),J.conj().T@model.apply(J))
        # Same parent coefficients also give full transfer at the common long
        # clock; here the axis columns are supplemented by the full K=4
        # all-input tests above, rather than mistaken for all references.
        budget.append(dict(label=label,eta=eta,curvature=eta*eta,
                           reference_exit_probability=probability,
                           finite_specified_rate=measured_rate,finite_worst_rate=worst_rate,
                           taylor_rate_allowance=remainder,
                           guaranteed_rate_interval=[eta-remainder,eta+remainder],
                           passes_curvature_budget=bool(eta<=beta),
                           axis_dimension=len(model.basis)))

    return dict(round=1325,status="PASS",new_axioms=0,count=len(checks),checks=checks,
                max_equality_residual=max(x["residual"] for x in checks if x["kind"]=="equality"),
                finite_role_examples=role_examples,
                small_complex_parent=small,full_unknown_transfers=transfer,
                coherent_witness=coherent,
                reference_budget=dict(L=L,nu=nu,beta=beta,tau=tau,
                                      common_full_transfer_time=float(common_time),
                                      source_contract="all coefficients fixed; only explicitly declared relative M_U orientation differs",
                                      cases=budget),
                limitations=[
                    "single-fermion finite full parent, not all many-body gauge fields",
                    "K>=4 role conclusion uses the full unknown-reference equivariant contract",
                    "full left coverage and the curvature budget are optional candidate tasks, not adopted axioms",
                    "coherent preparation, internal phase and readout resources are inputs",
                    "Duhamel remains sufficient; necessity of the Gram bound is only for the declared initial-curvature task",
                    "three/four/five flavor witnesses do not select the observed family number",
                    "mode energies and preparation responses are not physical mass poles or Higgs functional sources",
                ])


def run():
    result=compute()
    paths=[HERE/"proof.md",HERE/"check.py",STAGE/"_admission/round1325_admission.md",
           STAGE/"1324/proof.md",STAGE/"1324/check.py",STAGE/"1324/results.json",
           STAGE/"1323/proof.md"]
    result["source_sha256"]={p.relative_to(ROOT).as_posix():digest(p) for p in paths}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    result=run()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        same_result(result,json.loads(target.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1325,mode="write-new" if args.write_new else "read-only",
                          passed=True,count=result["count"],
                          max_equality_residual=result["max_equality_residual"],
                          actual_kappa=result["coherent_witness"]["actual_kappa"],
                          leakage=[x["reference_exit_probability"] for x in result["reference_budget"]["cases"]])))


if __name__=="__main__":
    main()
