"""Independent mass-eigenstate CP-odd one-Higgs couplings for the two-source flavor texture."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


def check_texture(y,source_real,source_complex,qr,qd):
    x=qr*source_real+qd*source_complex
    ins=np.zeros_like(y)
    ins[1,2]=x
    ul,sm,vh=np.linalg.svd(y)
    order=np.argsort(sm);sm=sm[order]
    ul=ul[:,order];ur=vh.conj().T[:,order]
    np.testing.assert_allclose(ul.conj().T@y@ur,np.diag(sm),atol=2e-12)
    # The baseline radial Y of the ordinary linear mass terms is included.
    reduced=np.diag(ul.conj().T@(y+ins)@ur)/sm
    p=reduced.imag
    lam=sm**2
    h=y@y.conj().T
    v=ins@y.conj().T
    j=(h[0,1]*h[1,2]*h[2,0]).imag/h[0,0].real
    first=(x*y[1,2].conjugate()).imag
    second=(h[1,1].real+abs(y[2,2])**2)*first+qd*j
    np.testing.assert_allclose([sum(p),lam@p,(lam**2)@p],[0,first,second],atol=3e-12)
    np.testing.assert_allclose(np.trace(np.linalg.solve(y,ins)).imag,0,atol=1e-13)
    np.testing.assert_allclose([np.trace(v).imag,np.trace(h@v).imag],[first,second],atol=3e-12)
    pred=[]
    for i in range(3):
        others=[jidx for jidx in range(3) if jidx!=i]
        denom=np.prod(lam[i]-lam[others])
        coeff=lam[i]-abs(y[0,1])**2-abs(y[2,1])**2
        pred.append((qd*j+coeff*first)/denom)
    np.testing.assert_allclose(pred,p,atol=3e-12,rtol=2e-9)
    # Exact inverse except if the interference or oriented cycle vanishes.
    inv_D=(lam@((lam-h[1,1].real-abs(y[2,2])**2)*p))/j
    inv_R=inv_D-(lam@p)/(source_real*y[1,2].imag)
    np.testing.assert_allclose([inv_R,inv_D],[qr,qd],atol=5e-8)
    # A finite Higgs background preserves det(Y) along this insertion.
    for k in (-.3,.37):
        np.testing.assert_allclose(np.linalg.det(y+k*ins)/np.linalg.det(y),1,atol=4e-12)
    return {"q_R_D":[qr,qd],"masses_in_yb_units":sm.tolist(),
            "p_d_s_b":p.tolist(),"scalar_even_d_s_b":reduced.real.tolist(),
            "moment1":float(first),"moment2_subtracted":float(qd*j),"oriented_cycle":float(j)}


def sample():
    z=read("../1489/results.json")["reconstruction"]
    a,b,c,e,f=z["normalized_real_edges_A_B_C_E_F"]
    real,coeff=z["two_real_coefficients"]
    part=coeff*np.exp(1j*math.radians(z["eta_example_deg"]))
    y=np.array([[0,a,0],[b,c,real+part],[0,e,f]],dtype=complex)
    np.testing.assert_allclose(y[1,2],complex(*z["Z_real_imag"]),atol=1e-13)
    return y,real,part


def arbitrary_edge_slopes(y,real,part):
    # A,B,C,E,F,R,D. The R and D columns represent parallel coherent paths.
    edges=[(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]
    values=[y[i,j] for i,j in edges[:5]]+[real,part]
    port=np.zeros((7,6))
    for edge,(i,j) in enumerate(edges):
        port[edge,i]=1
        port[edge,j+3]=1
    cycles=np.array([[0,0,0,0,0,-1,1],[0,0,-1,1,-1,0,1]])
    assert np.linalg.matrix_rank(port)==5 and np.linalg.matrix_rank(cycles)==2
    np.testing.assert_array_equal(cycles@port,np.zeros((2,6)))
    UL,sm,VH=np.linalg.svd(y);UR=VH.conj().T
    h=y@y.conj().T
    j=(h[0,1]*h[1,2]*h[2,0]).imag/h[0,0].real
    weights=np.array([real*y[1,2].imag,j])
    def direct(slopes):
        z=np.zeros_like(y)
        for idx,((i,k),val) in enumerate(zip(edges,values)):
            z[i,k]+=slopes[idx]*val
        coupling=np.diag(UL.conj().T@z@UR)/sm
        ps=coupling.imag
        l=sm**2
        obs=[l@ps,((l-h[1,1].real-abs(y[2,2])**2)*l)@ps]
        np.testing.assert_allclose(sum(ps),0,atol=2e-11)
        np.testing.assert_allclose(obs,weights*(cycles@slopes),atol=3e-12)
        return coupling
    response=np.array([direct(row).imag for row in np.eye(7)]).T
    assert np.linalg.matrix_rank(response,tol=1e-8)==2
    np.testing.assert_allclose(response@port,0,atol=1e-10)
    rng=np.random.default_rng(16481)
    for _ in range(25):direct(rng.normal(size=7))
    shift=direct(port@np.array([.3,-.8,1.1,.2,-.4,.5]))
    np.testing.assert_allclose(shift.imag,0,atol=1e-11)
    assert np.max(abs(shift.real))>.1
    return {"order":["A","B","C","E","F","R","D"],"two_cycles":cycles.tolist(),
            "port_rank":5,"p_rank":2,"port_example_real_scalar_shifts_heavy_to_light":shift.real.tolist()}


def finite_tests():
    rng=np.random.default_rng(1648)
    for _ in range(36):
        edges=rng.uniform(.22,1.5,5)*rng.choice([-1.,1.],5)
        a,b,c,e,f=edges
        real,mag=rng.uniform(.25,1.1,2)
        part=mag*np.exp(1j*rng.uniform(.4,2.6))
        y=np.array([[0,a,0],[b,c,real+part],[0,e,f]])
        qr,qd=rng.uniform(-2.1,2.1,2)
        r=check_texture(y,real,part,qr,qd)
        co=check_texture(y.conj(),real,part.conjugate(),qr,qd)
        np.testing.assert_allclose(np.array(co["p_d_s_b"]),-np.array(r["p_d_s_b"]),atol=1e-11)
    return 36


def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):return [stable(v) for v in x]
    if isinstance(x,float):return 0. if abs(x)<1e-14 else float(f"{x:.11g}")
    return x


def main():
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true")
    args=p.parse_args()
    y,r,d=sample()
    out={"independent_random_cases_and_conjugates":finite_tests(),
         "calibrated_fixture":[check_texture(y,r,d,qr,qd) for qr,qd in ((0,2),(2,0),(2,2),(0,0))],
         "all_seven_edges":arbitrary_edge_slopes(y,r,d)}
    eq=out["calibrated_fixture"][2]
    assert abs(eq["moment1"])<1e-13 and np.max(abs(np.array(eq["p_d_s_b"])))>1e-4
    out=stable(out)
    if args.check:assert out==read("independent_results.json")
    else:(HERE/"independent_results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for v in out["calibrated_fixture"]:print(v["q_R_D"],v["p_d_s_b"])

if __name__=="__main__":
    main()
