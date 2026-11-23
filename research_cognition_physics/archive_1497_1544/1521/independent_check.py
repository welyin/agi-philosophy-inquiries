"""Independent crossed scalar-current test in the complete 10 flavor space."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))

def flavor_frames():
    spec=importlib.util.spec_from_file_location("round1516",HERE.parent/"1516"/"check.py")
    src=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(src)
    core=src.parent_core(read("../1516/inputs.json"))
    spectra, down, lepton, hdabs, _, _ = core
    down_frame, roots, _=src.mixing_core(core,12)
    phase=np.exp(1j*math.radians(12))
    bare=hdabs[1:]/abs(1+2*phase*roots)
    he=np.r_[hdabs[0],bare*(1-3*phase*roots)]
    hd=np.r_[hdabs[0],bare*(1+2*phase*roots)]
    def gram(h):
        y=np.column_stack([np.zeros(3),np.diag(spectra)])
        return y@(np.eye(4)-np.outer(h,h.conj())/np.vdot(h,h))@y.T
    e_values, e_frame=np.linalg.eigh(gram(he))
    np.testing.assert_allclose(e_values,lepton**2,rtol=1e-8)
    np.testing.assert_allclose(down_frame.conj().T@gram(hd)@down_frame,np.diag(down**2),atol=3e-19)
    return down_frame,e_frame

def check():
    down,el=flavor_frames()
    up_spectrum=np.array(read("../1397/inputs.json")["high_scale"]["up"])[:,0]
    mass=np.array(read("../1518/results.json")["witness"]["even_triplets_GeV"])
    fits=read("../1516/results.json")["solutions"]
    # The unique unit ray common to the two light planes annihilates b and tau.
    _,sing,vh=np.linalg.svd(np.stack([down[:,2].conj(),el[:,2].conj()]))
    assert min(sing)>0.1
    ray=vh[-1].conj()
    dcoef=down[:,:2].conj().T@ray
    ecoef=el[:,:2].conj().T@ray
    np.testing.assert_allclose(down[:,:2]@dcoef,ray,atol=1e-13)
    np.testing.assert_allclose(el[:,:2]@ecoef,ray,atol=1e-13)
    l1=float(np.linalg.norm(dcoef,1)*np.linalg.norm(ecoef,1))
    rng=np.random.default_rng(1521)
    out={}
    for key,sol in fits.items():
        o=np.array(sol["12"]["one_passing_frame"]["O"])
        u=o[0]
        y=o.T@np.diag(up_spectrum)@o
        # Check once with actual complex independent heavy Yukawas, once with a
        # collinear up witness that makes cancellation of the holomorphic term possible.
        overlap=abs(np.vdot(u,ray))**2
        quad=up_spectrum[0]**2*overlap/mass[0]**2
        assert np.isclose(abs(u@y@ray)**2/(mass[0]**2*quad),1)
        for collinear in (True,False):
            blocks=[]
            for _ in range(2):
                t=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
                blocks.append((t+t.T)/2)
            up_all=[y, (mass[1]/mass[0])*y, np.zeros_like(y)] if collinear else [y]+blocks
            port1=np.stack([u@ya@down[:,:2]/m for ya,m in zip(up_all,mass)])
            port2=np.stack([u@ya@el[:,:2]/m for ya,m in zip(up_all,mass)])
            crossed=port1.T@port2.conj()
            shared=dcoef.T@crossed@ecoef.conj()
            projected=sum(abs(u@ya@ray/m)**2 for ya,m in zip(up_all,mass))
            np.testing.assert_allclose(shared,projected,rtol=1e-9,atol=1e-39)
            assert max(abs(crossed).flat)>quad/l1
            all1=np.stack([(o@ya@down).ravel()/m for ya,m in zip(up_all,mass)])
            all2=np.stack([(o@ya@el).ravel()/m for ya,m in zip(up_all,mass)])
            kk=np.stack([ya.ravel()/m for ya,m in zip(up_all,mass)])
            np.testing.assert_allclose(np.linalg.norm(all1.T@all2.conj())/np.linalg.norm(kk.T@kk.conj()),1)
            if collinear:
                same=np.array([1.,1.,0])
                other=np.array([1.,-1.,0])
                assert abs(sum(s*b for s,b in zip(same,other)))<1e-15
                collinear_max=float(abs(crossed).max())
        out[key]=dict(gamma_mod180=sol["12"]["one_passing_frame"]["gamma_deg"],
            shared_ray_projection_on_up_squared=float(overlap),
            four_physical_components_max_lower_GeVminus2=quad/l1,
            two_collinear_arms_example_max_GeVminus2=collinear_max)
    # A different right charged-lepton orientation in a wider class can hide u in tau.
    alt=np.eye(3)[:,[1,2,0]]
    assert np.linalg.norm(np.eye(3)[0]@np.diag(up_spectrum)@alt[:,:2])==0
    sin_beta_squared = math.sin(math.atan(read("inputs.json")["tan_beta"]))**2
    primary = read("results.json")["witnesses"]
    full, hidden = [primary[k] for k in ("shared_tensor_spends_full_trace","extra_tensor_hidden_from_up_quark")]
    def coherent_down_gain(branch):
        return 1+(branch["c_heavy_to_light"]*mass[0]/mass[1])**2
    rate_ratio = (coherent_down_gain(full)/coherent_down_gain(hidden))**2
    for tag in ("Gamma_pi_over_limit","Gamma_K_over_limit"):
        got = full["scalar_exchange_only_central"][tag]/hidden["scalar_exchange_only_central"][tag]
        assert abs(got/rate_ratio-1)<1e-10
    base_full_flavor = np.dot(up_spectrum,up_spectrum)/(9*mass[0]**2)
    for entry in out.values():
        entry["four_components_lower_at_benchmark_tan_beta"] = entry["four_physical_components_max_lower_GeVminus2"]/sin_beta_squared
    return dict(complete_81_max_reduced_coefficient_lower_GeVminus2=base_full_flavor,
        complete_81_lower_at_benchmark_tan_beta=base_full_flavor/sin_beta_squared,
        ratio_of_the_two_down_source_witness_partial_rates=rate_ratio,
        projective_l1_gain=l1, intersection_normal_singular_values=sing.tolist(),
        light_flavor_cases=out)

def trees(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            trees(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=0)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    m=p.add_mutually_exclusive_group(required=True)
    m.add_argument("--write",action="store_true")
    m.add_argument("--check",action="store_true")
    args=p.parse_args()
    values=check()
    if args.write:
        with (HERE/"independent_results.json").open("w",encoding="utf8") as h:
            json.dump(values,h,indent=2,allow_nan=False)
            h.write("\n")
    else:
        trees(values,read("independent_results.json"))
        print("1521 independent checks passed")
