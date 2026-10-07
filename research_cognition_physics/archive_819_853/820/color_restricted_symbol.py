"""820 working probe: color-momentum-only compensation at fixed geometry.
Uses the actual 753 constant color connection and electric background.
Fourier sampling is not an all-frequency smooth right-inverse proof.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
with ResearchRuntime(Layout()).installed():
    import joint_irreducible_source_completion as previous
bg=previous.bg
TARGET=HERE/'color_restricted_symbol_results.json'

def data():
    c,C=previous.connection_matrices()
    A=np.array([previous.components(a) for a in c['A']])
    E=np.array([previous.components(e) for e in c['E']])
    F=np.array([[previous.components(f) for f in row] for row in c['F']])
    return A,E,F,C

def symbol(k):
    A,E,F,C=data();I=np.eye(8);k=np.asarray(k,float)
    G=np.hstack([1j*k[i]*I+C[i] for i in range(3)])
    M=F.reshape(3,24);energy=E.reshape(1,24)
    return np.vstack((G,M,energy))

def run():
    A,E,F,C=data();zero=symbol([0,0,0])
    left=np.zeros((12,3));left[:8]=-A.T;left[8:11]=np.eye(3)
    zero_null_error=float(np.max(abs(zero.conj().T@left)))
    eigen0=np.linalg.eigvalsh(zero@zero.conj().T)
    assert zero_null_error<1e-14 and sum(eigen0<1e-12)==3
    rows=[];smallest=float('inf');rank_losses=[]
    for x in range(-4,5):
        for y in range(-4,5):
            for z in range(-4,5):
                k=np.array([x,y,z])
                if not np.any(k):continue
                B=symbol(k);e=np.linalg.eigvalsh(B@B.conj().T);smallest=min(smallest,float(e[0]))
                if e[0]<1e-12:rank_losses.append(k.tolist())
    for direction in ([1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1],[1,2,-3]):
        for size in (1,4,16,64,256):
            k=size*np.array(direction)
            B=symbol(k)
            # Project out Gauss rows before measuring the remaining algebraic rank.
            G=B[:8];Q=np.eye(24)-G.conj().T@np.linalg.solve(G@G.conj().T,G)
            reduced=B[8:]@Q@B[8:].conj().T
            eigen=np.linalg.eigvalsh(reduced)
            rows.append(dict(k=k.tolist(),reduced_min=float(eigen[0]),
                             reduced_eigenvalues=eigen.tolist()))
    return dict(round=820,formal_round_completed=False,source='Original753 fixed color data.',
        zero_mode_eigenvalues=eigen0.tolist(),zero_mode_three_momentum_identities_error=zero_null_error,
        zero_mode_compatibility='mean(M_i) - A_i dot mean(G) = 0',
        sampled_nonzero_modes=728,sampled_min_Gram_eigenvalue=smallest,
        sampled_rank_losses=rank_losses,large_frequency_rows=rows,
        continuous_all_mode_inverse_proven=False,actual819_input_compensated=False,
        geometry_fixed_full_model_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
