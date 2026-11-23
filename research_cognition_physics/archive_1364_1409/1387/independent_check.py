"""Hold the visible quark row and both valence weak charges fixed, vary strange flavor."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rotated(n, axis, heavy, loss):
    a, b = np.zeros(n), np.zeros(n)
    a[:len(axis)] = axis
    b[heavy] = 1.
    return (np.eye(n)+(sqrt(1-loss)-1)*(np.outer(a,a)+np.outer(b,b))
            +sqrt(loss)*(np.outer(b,a)-np.outer(a,b)))


def currents(row0, x, rotations):
    e = np.array([1.,0.,0.])
    v = e-row0
    v /= np.linalg.norm(v)
    n = 3+len(rotations)
    t = np.zeros((4,n))
    t[:3,:3] = np.eye(3)-2*np.outer(v,v)
    d = np.eye(n)
    for k,(axis,loss) in enumerate(rotations):
        d = d @ rotated(n,axis,3+k,loss)
    u = rotated(4,e,3,x)
    w = u.T @ t @ d
    up = w @ w.T
    down = w.T @ w
    q_u, q_d = up[0,0]-1, 1-down[0,0]
    return w[0,:3], np.array([2*q_u+q_d,188*q_u+211*q_d,abs(down[0,1])])


def calculate():
    history = json.loads((HERE.parent/"1385"/"independent_results.json").read_text(encoding="utf-8"))
    m = np.array(history["row_dsb"])
    d = history["loss"]
    v = m / np.linalg.norm(m)
    x = .0004
    y = (d-x)/(1-x)
    a = v[0]**2*y
    b = 1-m[1]**2/(1-x-m[0]**2/(1-a)-m[2]**2)
    assert 0 < a < 1 and 0 < b < 1
    first, charges1 = currents(v,x,[(v,y)])
    bare = m / np.sqrt([1-a,1-b,1.]) / sqrt(1-x)
    assert abs(np.linalg.norm(bare)-1)<1e-14
    second, charges2 = currents(bare,x,[(np.eye(3)[0],a),(np.eye(3)[1],b)])
    assert np.allclose([first,second],m[None,:],atol=1e-13,rtol=0)
    assert np.allclose(charges1[:2],charges2[:2],atol=1e-13,rtol=0)
    assert abs(charges1[0]-.0011) < np.linalg.norm([.0045,.0003])
    assert abs(charges1[1]-.65) < np.linalg.norm([.29,.32,.02])
    assert abs(charges2[2]) < 1e-13
    rate_hi = (9.6+2*1.9)*1e-11
    budget = (8.4+2.)*1e-11
    k_low = .04965-2*(.00038+.00037)
    vus_hi = (.2157+2*.0006)/(.9704-2*.0032)
    conversion = 3*.8*k_low/(2*vus_hi**2)
    cap = (sqrt(rate_hi)+sqrt(budget))/sqrt(conversion)
    assert charges1[2]>cap and 9.6e-11 < budget
    return {"up_singlet_probability":x, "common_diagonal_d_loss":a,
            "two_singlet_diagonal_s_loss":b, "one_singlet_projection":float(charges1[2]),
            "two_singlet_projection":float(charges2[2]),
            "matched_neutral_shifts":charges1[:2].tolist(), "rare_modulus_cap":cap,
            "minimum_other_rate_for_one":(sqrt(conversion)*charges1[2]-sqrt(rate_hi))**2}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    option=ap.parse_args()
    result=calculate()
    if option.check:
        old=json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        assert old.keys()==result.keys()
        for name in result:
            np.testing.assert_allclose(result[name],old[name],rtol=1e-10,atol=1e-15)
        print("1387: independent rate geometry and two explicit charged-neutral realizations verified.")
    else:
        print(json.dumps(result,indent=2))
