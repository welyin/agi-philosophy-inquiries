"""Born cross section for the registered common electroweak current and stable pairs."""
import math

import numpy as np


def window(masses, v, mh, gut_g1, g2):
    energy = 2000.0
    recoil_edges = (1100.0, 1300.0, 1500.0)
    gy = math.sqrt(3/5)*gut_g1
    mz2 = v**2*(g2**2 + gy**2)/4
    sw2 = gy**2/(gy**2 + g2**2)
    gf = 1/(math.sqrt(2)*v**2)
    cm2 = energy**2
    mass2 = mh**2
    unit_fb = 3.89379372e11
    x, weights = np.polynomial.legendre.leggauss(96)

    def produce(q):
        phase = (1-(np.sqrt(q)+math.sqrt(mz2))**2/cm2)
        phase *= 1-(np.sqrt(q)-math.sqrt(mz2))**2/cm2
        assert min(phase) > 0
        front = gf**2*mz2**2*(1+(4*sw2-1)**2)/(96*math.pi*cm2)
        return (unit_fb*front*np.sqrt(phase)*(phase+12*mz2/cm2)
                / (1-mz2/cm2)**2)

    rows = []
    for low, high in zip(recoil_edges, recoil_edges[1:]):
        half = (high**2-low**2)/2
        q = (high**2+low**2)/2 + x*half
        absorptive = sum(m**4*np.sqrt(np.maximum(0, 1-4*m**2/q))
                         for m in masses)/(8*math.pi*v**2)
        stripped = absorptive/(math.pi*(q-mass2)**2)
        rows.append({"window_GeV": [low, high],
                     "normalized_integral": float(half*(weights @ stripped)),
                     "rate_fb": float(half*(weights @ (stripped*produce(q))))})
    return {"sqrt_s_GeV": energy, "bins": rows,
            "normalized_integral": sum(row["normalized_integral"] for row in rows),
            "rate_fb": sum(row["rate_fb"] for row in rows)}
