"""Rational majorant retaining original group characters; full sphere moments."""
from fractions import Fraction as F
import itertools
import math
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import block_majorant_probe as blocks


def polynomials(n=16):
    poch=[1]
    for j in range(n):poch.append(poch[-1]*(5+j))
    out=[]
    for k in range(n//2+1):
        coeff=[F(0)]*(n-2*k+1)
        for ell in range(n//2-k+1):
            d=k+ell;remaining=n-2*d
            scalar=F(math.factorial(n)*4**k*(ell+1),2**n*20**ell)
            for r in range(remaining+1):
                s=remaining-r
                factor=scalar/(poch[r+d]*poch[s+d])
                f=[F(math.comb(k+2+i,i)*math.comb(ell+1+r-i,r-i))*F(2)**(2*i-r) for i in range(r+1)]
                g=[F(math.comb(k+2+j,j)*math.comb(ell+1+s-j,s-j))*F(2)**(2*j-s) for j in range(s+1)]
                for i,x in enumerate(f):
                    for j,y in enumerate(g):coeff[i+j]+=factor*x*y
        assert all(x>0 for x in coeff)
        out.append(coeff)
    return out


def evaluate(a,cs,polys):
    h=[np.ones_like(a)]+[np.zeros_like(a) for _ in range(8)]
    for c in cs:
        for k in range(1,9):h[k]=h[k]+c*h[k-1]
    out=np.zeros_like(a)
    for k,coef in enumerate(polys):
        v=np.zeros_like(a)
        for x in reversed(coef):v=v*a+float(x)
        out+=v*h[k]
    return out


def weights():
    result=[]
    for n in (0,2,4):
        for subset in itertools.combinations(range(5),n):
            w=tuple(int(j in subset)-int(4 in subset) for j in range(4))
            odd=sum(j in subset for j in (3,4))%2
            result.append((w,odd))
    return result


def grid_float(n=28):
    coords=np.array(np.unravel_index(np.arange(n**4),(n,)*4)).T
    roots=np.exp(2j*np.pi*np.arange(n)/n)
    rs=[(roots[(coords@np.array(w))%n],odd) for w,odd in weights()]
    a=(464-sum(2*r.real for r,odd in rs if odd))/1280
    cs=[(42+30*roots[coords[:,j]].real)/3200 for j in range(3)]
    value=evaluate(a,cs,polynomials())
    even=np.ones(n**4)
    for r,odd in rs:
        if not odd:even*=(2+2*r.real)/16
        else:even*=(19-r.real)/20
    z=[roots[coords[:,j]] for j in range(4)]+[roots[(-coords.sum(axis=1))%n]]
    density=np.ones(n**4)/12
    for i,j in ((0,1),(0,2),(1,2),(3,4)):density*=abs(z[i]-z[j])**2
    return float(np.mean(value*even*density)),float(value.max()),float(even.mean())


if __name__=='__main__':
    print('polynomials',sum(len(p) for p in polynomials()),flush=True)
    print('integral',grid_float(),flush=True)
