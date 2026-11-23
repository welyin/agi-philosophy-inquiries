"""Exact residue certificate for the low-energy flavor-observable Jacobian."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


class DualMatrix:
    """Matrix and directional derivative over one small prime field."""
    def __init__(self, value, tangent, p):
        self.p = p
        self.value = value % p
        self.tangent = tangent % p

    def __add__(self, other):
        return DualMatrix(self.value+other.value,self.tangent+other.tangent,self.p)

    def __matmul__(self, other):
        return DualMatrix(self.value@other.value,
                          self.tangent@other.value+self.value@other.tangent,self.p)

    def scale(self, x, dx=0):
        return DualMatrix(x*self.value,x*self.tangent+dx*self.value,self.p)

    def transpose(self):
        return DualMatrix(self.value.T,self.tangent.T,self.p)

    def centered(self):
        eye = np.eye(3,dtype=np.int64)
        factor = pow(3,-1,self.p)
        return DualMatrix(self.value-eye*(int(self.value.trace())*factor),
                          self.tangent-eye*(int(self.tangent.trace())*factor),self.p)

    def inverse(self):
        block = np.column_stack((self.value,np.eye(3,dtype=np.int64)))
        for i in range(3):
            pivot = next(k for k in range(i,3) if block[k,i])
            block[[i,pivot]] = block[[pivot,i]]
            block[i] = block[i]*pow(int(block[i,i]),-1,self.p)%self.p
            for k in range(3):
                if k!=i:
                    block[k] = (block[k]-block[k,i]*block[i])%self.p
        inv = block[:,3:]
        return DualMatrix(inv,-(inv@self.tangent%self.p)@inv,self.p)


def mass_forms(x, dx, imag, p, pairs):
    h = DualMatrix(np.diag(x[:3]),np.diag(dx[:3]),p)
    f, df = np.zeros((2,3,3),dtype=np.int64)
    for k,(i,j) in enumerate(pairs):
        f[i,j] = f[j,i] = x[3+k]+imag*x[9+k]
        df[i,j] = df[j,i] = dx[3+k]+imag*dx[9+k]
    f = DualMatrix(f,df,p)
    sf = f.scale(x[16]+imag*x[17],dx[16]+imag*dx[17])
    n = h+sf.scale(-3)
    return [(h+sf).scale(x[15],dx[15]),h+f,h+f.scale(-3),
            (n@f.inverse()@n).scale(x[18],dx[18])]


def response(x, dx, imag, p, pairs):
    minus = mass_forms(x,dx,-imag,p,pairs)
    plus = mass_forms(x,dx,imag,p,pairs)
    up,down,lepton,neutral = [a.transpose()@b for a,b in zip(minus,plus)]
    outputs = []
    for mass in (up,down,lepton):
        outputs.extend((mass,mass@mass,mass@mass@mass))
    def mixed(a,b):
        return [a@b,a@b@b,a@a@b,a@a@b@b]
    outputs.extend(mixed(up,down))
    centered = neutral.centered()
    outputs.extend((centered@centered,centered@centered@centered))
    outputs.extend(mixed(lepton,centered))
    return np.array([obj.tangent.trace()%p for obj in outputs],dtype=np.int64)


def determinant(matrix, p):
    work = matrix.copy()%p
    det = 1
    for i in range(len(work)):
        rows = np.flatnonzero(work[i:,i])
        if not len(rows):
            return 0
        j = i+int(rows[0])
        if i!=j:
            work[[i,j]] = work[[j,i]]
            det = -det
        pivot = int(work[i,i])
        det = det*pivot%p
        inverse = pow(pivot,-1,p)
        for k in range(i+1,len(work)):
            work[k] = (work[k]-(int(work[k,i])*inverse%p)*work[i])%p
    return int(det)


def verify_finite_graph(cfg):
    """Check the normalized light subspace before rounding its tiny masses."""
    f = np.zeros((3,3),dtype=complex)
    for k, (i,j) in enumerate(cfg["pairs"]):
        f[i,j] = f[j,i] = cfg["F_real_upper"][k]+1j*cfg["F_imag_upper"][k]
    h = np.diag(cfg["H_diag"])
    r, s, t = cfg["r"], complex(cfg["s_real"],cfg["s_imag"]), cfg["t"]
    direct = r*(h-3*s*f)
    right = -r*r*f/t
    leading = -direct@np.linalg.solve(right,direct.T)
    relative_errors = []
    for scale in (2.**-13,2.**-14):
        graph = -scale*np.linalg.solve(right,direct.T)
        for _ in range(8):
            graph = scale*np.linalg.solve(right,graph.conj()@direct@graph-direct.T)
        w, vec = np.linalg.eigh(graph.conj().T@graph)
        norm = (vec*(1/np.sqrt(1+w)))@vec.conj().T
        sub = np.vstack((norm,graph@norm))
        light = np.linalg.solve(norm.conj(),direct@graph@norm)
        full = np.block([[np.zeros((3,3)),direct],[direct.T,right/scale]])
        np.testing.assert_allclose(sub.conj().T@sub,np.eye(3),atol=1e-14)
        np.testing.assert_allclose(full@sub,sub.conj()@light,atol=1e-12)
        np.testing.assert_allclose(light,light.T,atol=1e-12)
        relative_errors.append(np.linalg.norm(light/scale-leading)/np.linalg.norm(leading))
    assert relative_errors[1] > 1e-11
    assert 3.95 < relative_errors[0]/relative_errors[1] < 4.05


def evaluate():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    verify_finite_graph(cfg)
    x = np.array(cfg["H_diag"]+cfg["F_real_upper"]+cfg["F_imag_upper"]+
                 [cfg["r"],cfg["s_real"],cfg["s_imag"],cfg["t"]],dtype=np.int64)
    receipts = {}
    for p in cfg["primes"]:
        assert p>3 and all(p%d for d in range(2,int(p**0.5)+1))
        imag = next(z for z in range(1,p) if z*z%p==p-1)
        basis = np.eye(len(x),dtype=np.int64)
        jac = np.column_stack([response(x,d,imag,p,cfg["pairs"]) for d in basis])
        det = determinant(jac,p)
        assert det
        # Reversing the representation of i conjugates all four tensors.
        flipped = np.column_stack([response(x,d,-imag,p,cfg["pairs"]) for d in basis])
        np.testing.assert_array_equal(jac,flipped)
        receipts[str(p)] = {"sqrt_minus_one":imag,"determinant_mod_p":det,"rank":19}
    return {"rational_witness_output_coordinates":len(x),"field_receipts":receipts,
            "physical_flavor_coordinates_including_majorana":22,
            "generic_local_codimension":3}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flag = parser.add_mutually_exclusive_group(required=True)
    flag.add_argument("--check",action="store_true")
    flag.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    out = evaluate()
    name = HERE/"results.json"
    if args.write_new:
        with name.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        assert json.loads(name.read_text(encoding="utf-8"))==out
        print("1554: nineteen independent real observables certified over both prime fields")
