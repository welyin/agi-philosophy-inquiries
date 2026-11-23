"""Independent shape and neutron axial checks at the archived 1757 charged-current points."""
import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent


def read(path):
    return json.loads(path.read_bytes())


def compare(point, rates):
    s, u = point["us"], point["ub"]
    d = math.sqrt(1 - s*s - u*u)
    dvec, sd = rates["holdout_nuclear"]
    svec, *errors = rates["holdout_Kl3"]
    ss = math.hypot(*errors)
    axial, sa = rates["holdout_axial"]

    # Positive real branch, fitting only the two vector normalizations.
    rd, rs = dvec - d, svec - s
    ad, a_s = 2*d - dvec, 2*s - svec
    axial_pred = a_s / ad
    neutron_factor = ad / dvec
    r = neutron_factor * svec / a_s
    assert math.isclose(r, svec / (dvec * axial_pred))
    # Same vector data, arbitrary complex phases: axial norms are distances
    # from the moving vector endpoint to the fixed point twice the left current.
    # Hold the down vector real, choose the positive imaginary sign for strange.
    cosine = (4*s*s + svec*svec - (axial*ad)**2) / (4*s*svec)
    assert -1 < cosine < 1
    vs = svec * complex(cosine, math.sqrt(1 - cosine**2))
    strange_right = vs - s
    assert math.isclose(abs(s - strange_right) / ad, axial)

    def r_interval(n):
        return [(svec-n*ss)/((dvec+n*sd)*(axial+n*sa)),
                (svec+n*ss)/((dvec-n*sd)*(axial-n*sa))]

    # Eq (20) in KTeV 0912.1291, and the CT discrepancy in Eq (12).
    rate_combination, width_B = 1.2407, .0044
    delta, width_delta = -.0035, .0080
    log_measure, width_log = .1915, .0122
    ratio_from_all_rates = svec / (dvec * axial)
    def shape(value):
        return math.log(rate_combination*value + delta)
    def log_box(n, low, high):
        return [math.log((rate_combination-n*width_B)*low+delta-n*width_delta),
                math.log((rate_combination+n*width_B)*high+delta+n*width_delta)]
    one = r_interval(1)
    # If the *same* rate products and the same hadronic references are
    # used in both factors, their product is exactly the hadronic ratio.
    roundtrip = []
    h, f = 1.195, .97  # illustrative normalization, cancels algebraically
    for point_axial, point_svec, point_dvec in [(axial,svec,dvec),(.235,.219,.976)]:
        observed_B = (h*point_axial)*point_dvec/(f*point_svec)
        reconstructed_r = point_svec/(point_dvec*point_axial)
        roundtrip.append(observed_B*reconstructed_r)
        assert math.isclose(roundtrip[-1],h/f)
    # The quoted Perkeo value is extracted under the real V,A hypothesis.
    # Inclusive amplitude convention, with the 2022 chiral estimate shown
    # as a range rather than a measured nuisance distribution.
    perkeo = 1.27641
    correction = [.014,.026]
    g_needed = [perkeo/((1+t)*neutron_factor) for t in reversed(correction)]
    g_SM = [perkeo/(1+t) for t in reversed(correction)]
    return dict(left=dict(d=d,s=s,u=u),positive_real_vector_fit=dict(Rd=rd,Rs=rs,
        withheld_axial_ratio=axial_pred,offset_in_axial_widths=(axial_pred-axial)/sa,
        r=r,predicted_logC=shape(r),neutron_VA_factor=neutron_factor,
        neutron_gA_QCD_required=g_needed,pure_left_gA_QCD_required=g_SM,
        fraction_shift_inferred_gA=1/neutron_factor-1),
        arbitrary_phase_central_example=dict(Rs_real=strange_right.real,
        Rs_imag=strange_right.imag,Rd=rd,axial=abs(s-strange_right)/ad),
        calibrated_all_three=dict(r=ratio_from_all_rates,r_one_width_box=one,
        logC=shape(ratio_from_all_rates),logC_one_width_box=log_box(1,*one),
        pure_left_logC=shape(1),pure_left_one_width_box=log_box(1,1,1),
        KTeV_measured_one_width=[log_measure-width_log,log_measure+width_log],
        same_input_product_witness=roundtrip,hadronic_only_reference=h/f))


def calculate():
    frozen = read(ARCH / "1757/results.json")["central"]
    independent = read(ARCH / "1757/independent_results.json")["combined_central_solution"]
    points = dict(fixed_source=frozen,independent_clocks={"us":independent["Vus"],"ub":independent["Vub"]})
    rates = read(ARCH / "1755/inputs.json")
    return {name:compare(p,rates) for name,p in points.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = calculate()
    output = HERE / "independent_results.json"
    if args.check:
        assert result == read(output)
        print("1759 CT shape, calibration cancellation and neutron axial matching verified")
    else:
        with output.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write("\n")
