"""725 entry: original wall reaches the radial axis; not a deficiency proof."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_material_boundary_transport as wall
import joint_quantum_reference_forms as old
TARGET=HERE/'wall_domain_entry_results.json'


def run():
    _,lam,_=wall.base_data();a=-lam;M=old.M;R2=6*M
    # Actual X_f for h>0 in the radial plane. Its smooth signed extension is
    # used only to bracket the first axis hit, not as a physical continuation.
    def vector(z):
        h,s=z;F=M-(h*h+s*s)/6;f=s+a*h
        return F*np.array([a-h*f/R2,1-s*f/R2])
    hmax=.05;smax=.1;h0=.04
    Fmin=M-(hmax*hmax+smax*smax)/6
    lower=Fmin*(a*(1-hmax*hmax/R2)-hmax*smax/R2)
    sbound=M*(1+smax*(smax+a*hmax)/R2)
    hitbound=h0/lower
    assert a>0 and lower>0 and sbound*hitbound<smax
    z=np.array([h0,0.]);dt=hitbound/1000;t=0.
    for _ in range(1100):
        previous=z.copy();before=t
        k1=-vector(z);k2=-vector(z+dt*k1/2)
        k3=-vector(z+dt*k2/2);k4=-vector(z+dt*k3)
        z=z+dt*(k1+2*k2+2*k3+k4)/6;t+=dt
        if z[0]<=0:break
    assert previous[0]>0>=z[0] and t<=hitbound+dt
    assert max(abs(previous[1]),abs(z[1]))<smax
    # Test the symmetrized chain rule on a local smooth radial function.
    def psi(h,s):return np.exp(-h*h-.3*s*s)
    def differential(which,h,s,u):
        F=M-(h*h+s*s)/6
        if which=='T':
            x=F*np.array([h*(1-h*h/R2),-s*h*h/R2])
            lap=F*(4-h*h/(2*M))
        elif which=='s':
            x=F*np.array([-h*s/R2,1-s*s/R2]);lap=-F*s/(3*M)
        else:
            x=vector(np.array([h,s]));lap=F*(-(s+a*h)/(3*M)+3*a/h)
        step=2e-6
        du=np.array([(u(h+step,s)-u(h-step,s))/(2*step),
                     (u(h,s+step)-u(h,s-step))/(2*step)])
        return -1j*old.HBAR*(x@du+lap*u(h,s)/2)
    errors=[]
    for h,s in ((.2,-.3),(.7,.5),(1.1,-.2)):
        direct=differential('f',h,s,psi)
        mapped=differential('s',h,s,psi)-lam*.5*(
            differential('T',h,s,lambda u,v:psi(u,v)/u)+differential('T',h,s,psi)/h)
        error=float(abs(direct-mapped)/(1+abs(direct)));errors.append(error)
        assert error<1e-8
    deps=('research_note_556.md','research_note_563.md','research_note_652.md','research_note_724.md',
          'joint_material_boundary_transport.py','joint_quantum_reference_forms.py')
    return dict(entry_round=725,new_formal_round=False,tests_run=2,failures=0,errors=0,
                result=dict(original_lambda=float(lam),radial_speed_lower=float(lower),
                    singlet_speed_upper=float(sbound),analytic_backward_hit_upper=float(hitbound),
                    numerical_hit_time_bracket=[float(before),float(t)],
                    signed_h_bracket=[float(previous[0]),float(z[0])],
                    s_at_crossing=float(z[1]),local_symmetrized_chain_rule_errors=errors,
                    original_complete_smooth_flow_argument_does_not_apply=True,
                    no_selfadjoint_extension_impossibility_claim=True),
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
                scope='Audit only: inherited radial-axis issue in the actual651 wall. Finite axis hit invalidates direct reuse of652 complete smooth flow proof, not all self-adjoint realizations or square-form observables. Next compare the closed first-order form with the same full energy and actual region records.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
