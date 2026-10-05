"""807 working: original principal symbol after Hadamard cross compression.

This does not construct a coupled on-shell boson or evaluate a switched signal.
It checks the previously missing frequency/sterile compression of a shear jet.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as original
TARGET=HERE/'sterile_principal_cross_probe_results.json'

def run():
    gx,gy,gz=original.original.GAMMA
    p=(np.eye(64)+gy)/2;q=np.eye(64)-p
    rows=[]
    for a in (1.,1j,.6+.8j):
        # 802 on-shell shear density corresponds to K principal = a*k_y*Gamma_z/2.
        k=.5*a*gz;cross=p@k@q
        response=(cross-cross.conj().T)/1j
        sterile=response[np.ix_([30,31],[30,31])]
        val=np.linalg.eigvalsh(sterile)
        assert abs(val[0]+abs(a)/2)<1e-13 and abs(val[1]-abs(a)/2)<1e-13
        rows.append(dict(shear=[float(np.real(a)),float(np.imag(a))],
            sterile_principal_response_eigenvalues=val.tolist()))
    conformal=p@gy@q
    assert np.max(abs(conformal))<1e-14
    return dict(working_round=807,all_checks_passed=True,original_Nambu_dimension=64,
        particle_sterile_indices=[30,31],three_complex_shear_checks=rows,
        conformal_cross_residual=float(np.max(abs(conformal))),
        Hadamard_leading_projector_only=True,boson_jet_proven_realizable=False,
        switched_continuous_response_evaluated=False,formal_round_completed=False,
        new_numbered_test_groups=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
