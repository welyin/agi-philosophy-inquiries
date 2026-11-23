"""Check the common mass columns through a Gram pencil and the heavy left projector."""
import argparse
import json
from math import log
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def data():
    a=np.sqrt([.0017,.0005])
    mass=np.array([.004,.075,3.,1500.,1800.])
    mu=.002
    eye=np.eye(2)

    def rows(amp):
        v=np.array([(.97373+amp)/2,(.22309+.23131*amp)/2])
        right=(np.array([.97373,.22309])-v)/a
        b=np.vstack((np.diag(right),np.zeros((1,2)),np.diag(np.sqrt(1-right**2))))
        column=mass[:,None]*b
        heavy=(mu*a - v@column[:2])/np.diag(column[-2:])
        first=np.array([*v,.0041*np.exp(-1.139j),*heavy])
        return first,b,column

    def residual(x):
        row=rows(x)[0]
        return np.vdot(row,row).real-1

    amp=.975
    for _ in range(6):
        amp-=residual(amp)/((residual(amp+1e-7)-residual(amp-1e-7))/2e-7)
    first,b,col=rows(amp)
    assert abs(residual(amp))<4e-15
    np.testing.assert_allclose(first@col,mu*a,atol=1e-15)
    clean=col-np.outer(first.conj(),first@col)
    q,_=np.linalg.qr(clean)
    leakage=np.linalg.norm(q[:2],axis=1)
    cc_sum_cap=float(np.prod(leakage))
    # A generalized Hermitian eigenproblem using a triangular Cholesky factor.
    metric=eye-np.outer(a,a)
    inverse=np.linalg.inv(np.linalg.cholesky(metric))
    pencil=inverse@(col.T@col-mu**2*np.outer(a,a))@inverse.T
    masses=np.sqrt(np.linalg.eigvalsh(pencil))
    right=a@b.T
    left_product=(first[0]*first[1].conj()).real
    assert abs(b[0]@b[1])<1.e-16
    assert abs(right[0]**2/(b[0]@b[0])+right[1]**2/(b[1]@b[1])-a@a)<1.e-16
    def std_loop(x):
        return x/8*((x+2)/(x-1)+3*(x-2)*log(x)/(x-1)**2)
    loop_cap=cc_sum_cap*max(std_loop((m/80.379)**2) for m in masses)/.225**5
    return {"t":amp,"heavy_mass_GeV":masses.tolist(),"same_bare_gram":(col.T@col).tolist(),
            "left_d_s_leak_norm":leakage.tolist(),"sum_heavy_abs_left_sd_upper":cc_sum_cap,
            "if_ordinary_X0_only_heavy_amplitude_cap":loop_cap,
            "min_sum_abs_heavy_right_sd":float(abs(right[0]*right[1])),
            "Fds":float(b[0]@b[1]),"left_light_u_lambda_re":left_product}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write-new",action="store_true")
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=data()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        original=json.loads(path.read_text(encoding="utf-8"))
        for k,value in result.items():
            np.testing.assert_allclose(value,original[k],atol=1e-18,rtol=2e-10)
        print("1465 independent OK")
