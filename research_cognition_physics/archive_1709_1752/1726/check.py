"""One common adjoint scale for the flavour insertion, X mass and radial decays."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
def read(file):
    return json.loads((HERE / file).read_bytes())


def compute():
    cfg = read("inputs.json")
    flavour = read(cfg["matter"])
    proton = read(cfg["proton"])
    cutoff = read(cfg["cutoff"])["cutoff_GeV"]
    v = read(cfg["electroweak_vev"])["v246_GeV"]
    m = cfg["illustrative_scalar_mass_GeV"]
    g = math.sqrt(proton["gauge_g_squared"])
    crossing = proton["gauge_crossing_GeV"]
    # Direct kinetic check for an off-diagonal real X generator with Tr(T*T)=1/2.
    vacuum = np.diag([2.,2.,2.,-3.,-3.])
    generator = np.zeros((5,5))
    generator[0,3] = generator[3,0] = .5
    comm = generator@vacuum - vacuum@generator
    radial_over_V = math.sqrt(2*np.trace(vacuum@vacuum))
    mx_over_gV = math.sqrt(2*np.trace(comm.T@comm))
    np.testing.assert_allclose([radial_over_V,mx_over_gV],[math.sqrt(60),5])
    f_at_crossing = crossing * radial_over_V / (mx_over_gV*g)
    rad = flavour["radial_flavour_test"]
    pair_norm_product = (rad["geometric_mean_of_dimensionless_offdiag_norms"]**2
                         * math.sqrt(rad["fraction_bs_among_all_flavor_changing_down"]
                                     *rad["fraction_tau_mu_among_all_flavor_changing_charged_leptons"]))
    # Real-scalar decay into two different Dirac families, summed over charge orders.
    rate_prefactor = math.sqrt(3)/(8*math.pi)*pair_norm_product*v**2/2
    at_crossing = rate_prefactor/f_at_crossing**2
    rows = {}
    for key, target in flavour["exact_and_conditional_implementations"].items():
        c = target["lift_singular_values"][0]
        canon = c*g*cutoff/crossing
        ratio = target["fraction_tagged_s_b_over_tau_mu"]
        rows[key] = dict(correction_norm=c,
                         Wilson_norm_at_reduced_Planck_and_crossing=canon,
                         cutoff_ceiling_at_crossing_for_Wilson_sqrt4pi=math.sqrt(4*math.pi)*crossing/(g*c),
                         MX_floor_at_reduced_Planck_for_Wilson_sqrt4pi=g*c*cutoff/math.sqrt(4*math.pi),
                         bs_and_tau_mu_widths_GeV_at_crossing_unit_radial_upper_vevs=
                            [m*at_crossing*math.sqrt(ratio),m*at_crossing/math.sqrt(ratio)])
    separate = {}
    f_amplitude = proton["amplitude_minimum_diagonal_phases"]
    for key in ("central_lattice", "one_width_each_lower_lattice"):
        w = proton[key]
        max_g_m_squared = w["low_energy_single_helicity_coefficient_upper_GeVm2"]/f_amplitude
        inv_f2_times_A = (mx_over_gV/radial_over_V)**2*max_g_m_squared
        np.testing.assert_allclose(inv_f2_times_A,
                                   (mx_over_gV*g/(radial_over_V*w["minimum_MX_GeV_at_unit_transfer"]))**2)
        separate[key] = dict(A_times_geometric_width_over_mass_upper=rate_prefactor*inv_f2_times_A,
                             A_times_geometric_width_GeV_at_example_mass_upper=m*rate_prefactor*inv_f2_times_A)
    return dict(real_radial_f_over_V=radial_over_V,MX_over_gV=mx_over_gV,
                canonical_f_GeV_at_crossing=f_at_crossing,product_of_pair_derivative_norms=pair_norm_product,
                geometric_width_over_mass_at_crossing_upper=at_crossing,
                geometric_width_GeV_at_example_mass_and_crossing_upper=m*at_crossing,
                matching_scale_GeV=cutoff,conditional_rows=rows,proton_to_flavour_bound=separate)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()), parse_float=lambda x: float(f"{float(x):.10g}"))
    path = HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
        print("1726 verified common adjoint normalization, rare-width product and ultraviolet coefficient costs.")
    else:
        with path.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
