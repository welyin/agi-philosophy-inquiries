"""Tag the charm row in the bottom-clock inverse problem with muonic D decays."""
import argparse
import json
import math
import runpy
from pathlib import Path


HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OLD = runpy.run_path(str(STAGE/'1766/check.py'))
CLOCK = OLD['CLOCKS']
PHASE = OLD['PHASE']
CS = OLD['BASE']['CFG']['mixing_coefficient_ps_inv']
D_MU = [48.02, .48, .24, .12, .15]       # fD |Vcd| in MeV, BESIII v2, corrected
DS_MU = [241.8, 2.5, 2.2]                # fDs |Vcs| in MeV, BESIII Born extraction
FD = [212.7, .6]                          # FNAL/MILC, MeV
FS_OVER_FD = [1.1749, .0016]              # FNAL/MILC, charged mass prescription
DS_RATE = [.993, .006]                    # long- plus short-distance width estimate


def band(data, n=1):
    width = n*math.hypot(*data[1:])
    return data[0]-width, data[0]+width


def charm(n):
    xd, xs = band(D_MU,n), band(DS_MU,n)
    fd, f_ratio, kappa = band(FD,n), band(FS_OVER_FD,n), band(DS_RATE,n)
    cd = (xd[0]/fd[1], xd[1]/fd[0])
    ratio = (xd[0]*f_ratio[0]*math.sqrt(kappa[0])/xs[1],
             xd[1]*f_ratio[1]*math.sqrt(kappa[1])/xs[0])
    return cd,ratio


def certificate():
    # Standard theory and physical labels exactly as in the 1766 bottom-clock inverse.
    md,ms = band(CLOCK['delta_md']),band(CLOCK['delta_ms'])
    bd,bs = band(CLOCK['mass_d']),band(CLOCK['mass_s'])
    xi_hi = max(band(CLOCK[label])[1] for label in ('FNAL_MILC','HPQCD'))
    q_hi = xi_hi*math.sqrt(md[1]/ms[0]*bs[1]/bd[0])
    pd2_hi = xi_hi**2*md[1]*bs[1]/(band(CS)[0]*bd[0])
    ps2_hi = ms[1]/band(CS)[0]
    h = pd2_hi+ps2_hi
    t_small_hi = 2*h/(1+math.sqrt(1-4*h))
    # d=|Vcd|. The angle ranges are in the first two quadrants;
    # all sin(beta+gamma) extrema below are taken on their actual arcs.
    d_lo,d_hi = charm(2)[0]
    ud2_lo = 1-d_hi**2-q_hi**2/(1+q_hi**2)
    gl,gh = map(math.radians,PHASE['68'])
    sm,sp = band(CLOCK['golden_S'])
    bl,bh = math.asin(sm)/2, math.asin(sp)/2
    rl,rh = math.pi/2-bh,math.pi/2-bl
    s_min = min(math.sin(bl+gl),math.sin(rh+gh))
    gsin = math.sin(gl)
    total_hi = pd2_hi/gsin**2*(1/ud2_lo+1/d_lo**2)
    ub2_over_cb2_hi = d_hi**2/(ud2_lo*s_min**2)
    assert 0 < bl+gl < math.pi/2 < rh+gh < math.pi
    assert total_hi < 1-t_small_hi and ub2_over_cb2_hi < 1
    assert math.pi/2 < rl+gl < rh+gh
    # Multiply the necessary column norm by t/(1+q^-2).
    # Using t<=1 gives an outward bound throughout the reflected-beta branch.
    reflected_hi = q_hi**2/(1+q_hi**2)/gsin**2 * (
        math.sin(rh)**2/ud2_lo + math.sin(rl+gl)**2/d_lo**2)
    assert reflected_hi < 1
    return dict(q_upper=q_hi,p_d_squared_upper=pd2_hi,small_t_upper=t_small_hi,
                Vud_squared_lower=ud2_lo, predicted_Vub2_plus_Vcb2_upper=total_hi,
                small_t_required_Vub2_plus_Vcb2_lower=1-t_small_hi,
                Vub2_over_Vcb2_upper=ub2_over_cb2_hi,
                reflected_beta_normalized_bottom_norm_upper=reflected_hi)


def cabibbo():
    # Wcb is an on-shell 13 TeV CMS determination with the experiment's template.
    cb_hi = .0459+2*math.hypot(.0044,.0064)
    vtb_lo = OLD['CFG']['top_decay']['Vtb_lower_95']
    rho_lo,rho_hi = charm(2)[1]
    previous = json.loads((STAGE/'1755/inputs.json').read_bytes())['inclusive_ub']
    out = {}
    for name,u in (('top_decay_only',math.sqrt(1-vtb_lo**2)),
                   ('with_inclusive_B_to_Xu',band(previous,2)[1])):
        c13 = math.sqrt(1-u*u)
        s23 = cb_hi/c13
        psi = math.atan(u*s23/math.sqrt(1-s23**2))
        lo = c13*math.sin(math.atan(rho_lo)-psi)
        hi = math.sin(math.atan(rho_hi)+psi)
        out[name] = dict(Vub_upper=u,Vcb_upper=cb_hi,psi_upper_rad=psi,
                         all_CKM_phases_Vus_outer_interval=[lo,hi])
    return out


def calculate():
    old = json.loads((STAGE/'1766/results.json').read_bytes())
    endpoints = json.loads((STAGE/'1766/independent_results.json').read_bytes())['refits']
    (d_lo,d_hi),(rho_lo,rho_hi) = charm(2)
    roots = {}
    for name, family in old['unitary_competitors'].items():
        xi = next(CLOCK[key][0] for key in ('HPQCD','FNAL_MILC') if key in name)
        pd = xi*math.sqrt(CLOCK['delta_md'][0]/CS[0]*CLOCK['mass_s'][0]/CLOCK['mass_d'][0])
        table = []
        for j,row in enumerate(family):
            be,ga = map(math.radians,(row['beta'],row['gamma']))
            w = row['Vcb']
            d = pd*math.sin(be+ga)/(math.sin(ga)*w)
            ratio = d/math.sqrt(1-w*w-d*d)
            table.append(dict(index=j,Vcd=d,cd_over_cs=ratio,Vcb=w,
                              charm_ratio_two_widths=rho_lo<=ratio<=rho_hi,
                              absolute_charm_two_widths=d_lo<=d<=d_hi))
        roots[name] = table
    real_material = {}
    for name,row in endpoints.items():
        s,u,w = (row['moduli'][key] for key in ('Vus','Vub','Vcb'))
        g,be = math.radians(row['gamma']),math.asin(row['S'])/2
        d = math.sqrt(1-u*u-s*s)*u*math.sin(be+g)/(w*math.sin(be))
        ratio = d/math.sqrt(1-d*d-w*w)
        real_material[name] = dict(Vcd=d,cd_over_cs=ratio,
                                   absolute_charm_two_widths=d_lo<=d<=d_hi,
                                   charm_ratio_two_widths=rho_lo<=ratio<=rho_hi)
    # Flavour-tagged inclusive W decays. Allow the heaviest final down quark
    # up to 5 GeV as a sensitivity case; neglect remaining masses here.
    r = (5./80.)**2
    zb = 1-(1-r)**2*(1+r/2)
    return dict(charm_readings={str(n):dict(Vcd=charm(n)[0],cd_over_cs=charm(n)[1])
                                 for n in (0,1,2)}, rigorous_box=certificate(),
                clock_roots=roots,six_charged_material_points=real_material,
                cabibbo_from_charm=cabibbo(),
                inclusive_RcW=dict(observed=[.489,.005,.019],massless=.5,
                                   maximal_tree_b_mass_shift=zb/(2*(2-zb))))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    checking=parser.parse_args().check
    results=json.loads(json.dumps(calculate()),parse_float=lambda x:float(f'{float(x):.10g}'))
    path=HERE/'independent_results.json'
    if checking:
        assert json.loads(path.read_bytes()) == results
        print('1772 independent OK: axial charm ratio, branch bounds, six material points and Cabibbo interval')
    else:
        with path.open('x',encoding='utf-8') as stream:
            json.dump(results,stream,indent=2,allow_nan=False)
            stream.write('\n')
