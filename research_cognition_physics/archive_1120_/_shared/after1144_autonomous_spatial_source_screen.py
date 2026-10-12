"""Short nonrelativistic reference screen. No parameter optimization/numbered round."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
FROZEN={'research_cognition_physics/archive_1120_/research_note_1128.md': '2c424638a38e87bd42234c3e6c3b77af860f2d81109be0dafa1e4c36c4f8cb44', 'research_cognition_physics/archive_1120_/research_note_1129.md': '2a3c198113b36ca17f7ebd19819ca9b0d14ceea7d5ad8b5b9c6ff3fad81d2fe0', 'research_cognition_physics/archive_1120_/research_note_1132.md': '4ebe32cc928cbae7365e7a69403bb3d7473d318983539fa701242d66cde6f624', 'research_cognition_physics/archive_1106_1119/research_note_1119.md': '931474ff6dfb64246269d00e075e0f0e2e6a4ff0fe6b8a2855e14fb3bb75e96c', 'research_cognition_physics/archive_1120_/research_note_1130.md': '4dec9307d0fc704c6ab174295d881f02b098f4d65c17e9ba1a2b0bdfc1a4d2a3', 'research_cognition_physics/archive_1120_/research_note_1135.md': '8e41e432e7a68b2244904c9291f8bc24936286185512e1d4133d1471787fdccb', 'research_cognition_physics/archive_1120_/1144/proof.md': '336ba9ab0cd448425e2704d3ec7ec6e683665369cace817ef520ad1f7ccd4ce2'}
OUT=Path(__file__).with_suffix('.json')
def main():
    for path,sha in FROZEN.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    T=10**9;sqrtTlower=31000
    assert sqrtTlower**2<T
    hcoef=2*F(110)*F(45,8)*F(7,4)/F(7,5)
    hbound=hcoef/sqrtTlower
    assert hbound<F(1,20)
    coherent_multiplier=F(5,7)+F(10,7)*F(1,20)
    field_bound=F(7,4)*coherent_multiplier
    grad_profile_squared=F(64,7)**2*F(88,21)
    gradient_bound=20*coherent_multiplier
    assert grad_profile_squared<400
    assert field_bound<2 and gradient_bound<16
    assert F(15,8)<2 and F(10,3)*F(7,4)<6
    gamma=F(5,4)
    boost_clock_rate_gap=1-1/gamma
    assert boost_clock_rate_gap==F(1,5)
    data={'status':'PASS_SHORT_REFERENCE_ONLY','numbered_round':None,'mathematical_checks':2,
      'new_physical_calibration':0,'frozen_sources':FROZEN,
      'normal_reference_constants':{'same_1144_T':T,'uniform_branch_displacement_upper':str(hbound),
        'coherent_multiplier_upper':str(coherent_multiplier),'translated_field_vector_upper':str(field_bound),
        'profile_gradient_L2_squared_upper':str(grad_profile_squared),'translated_field_gradient_vector_upper':str(gradient_bound),
        'clock_integrand_coefficient':'12*abs(lambda)/T^(3/2)',
        'position_integrand_coefficient':'32*sqrt(3)*abs(lambda)/sqrt(T)'},
      'old_bare_clock_connection':{'gamma':str(gamma),'proper_rate':'4/5','direct_sum_internal_rate':'1',
        'rate_gap_if_identified_with_same_proper_clock':str(boost_clock_rate_gap),
        'gap_removed_by_increasing_inertia':False,'independent_new_general_clock_no_go':False},
      'scope':{'autonomous_normal_reference_full_QR_bound_analytic':True,
        'same_total_energy_and_spatial_translation_analytic':True,
        'free_controller_and_position_parameters_optimized':False,
        'same_old_SR_clock_identity':False,'Lorentz_four_momentum':False,
        'point_material_microcausality':False,'complete_goal_candidate':False}}
    if '--write' in sys.argv:
        assert not OUT.exists(),'Do not overwrite saved screen'
        OUT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(OUT.read_text(encoding='utf-8'))
    print(json.dumps({'status':data['status'],'checks':2,'same_1144_field_displacement_upper':float(hbound),
        'no_new_numbered_round':True},ensure_ascii=False))
if __name__=='__main__':main()
