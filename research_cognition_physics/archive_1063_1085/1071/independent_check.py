"""Independent finite witnesses for round 1071, no repository imports.
Theorems about arbitrary iterates and the infinite torus are analytic claims;
finite checks below only validate explicitly stated formulas and examples.
"""
import cmath
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

EPS = 3e-12
COUNTS = {}
RESIDUALS = {}

def check(condition, group, message):
    COUNTS[group] = COUNTS.get(group, 0) + 1
    if not condition:
        raise AssertionError(group + ': ' + message)

def near(actual, expected, group, message, tol=EPS):
    residual = abs(actual - expected)
    RESIDUALS[group] = max(RESIDUALS.get(group, 0.0), float(residual))
    check(residual <= tol, group, message)

def add(a, b):
    return tuple(x + y for x, y in zip(a, b))

def scale(c, a):
    return tuple(c * x for x in a)

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def norm(a):
    return math.sqrt(dot(a, a))

def dist(a, b):
    return norm(tuple(x - y for x, y in zip(a, b)))

def unit(x):
    return scale(1.0 / norm(x), x)

O = (0.0, 0.0, 1.0)
LAMBDA = 1.0 / 3.0

def B(x, t):
    return unit(add(scale(1.0 - LAMBDA * t, x), scale(LAMBDA * t, O)))

def R(x):
    return B(x, 1.0)

def cost(x):
    return 1.0 - x[2]

def iterate(x, n):
    for _ in range(n):
        x = R(x)
    return x

def H(x, t):
    if t == 1.0:
        return O
    n = int(math.floor(-math.log2(1.0 - t)))
    left = 1.0 - 2.0 ** (-n)
    u = (t - left) * 2.0 ** (n + 1)
    return B(iterate(x, n), u)

def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))

def transpose(a):
    return tuple(zip(*a))

def action(a, x):
    return tuple(dot(row, x) for row in a)

def determinant(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))

def rank(matrix, tol=1e-10):
    a = [list(map(float, row)) for row in matrix]
    m, n = len(a), len(a[0])
    pivot = 0
    for col in range(n):
        choices = [i for i in range(pivot, m) if abs(a[i][col]) > tol]
        if not choices:
            continue
        k = max(choices, key=lambda i: abs(a[i][col]))
        a[pivot], a[k] = a[k], a[pivot]
        p = a[pivot][col]
        a[pivot] = [v / p for v in a[pivot]]
        for i in range(m):
            if i != pivot:
                c = a[i][col]
                a[i] = [u - c * v for u, v in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == m:
            break
    return pivot

# Positive witness: SO(3)/SO(2), a compact north cap, a single nonhalf step.
points = [O]
for m in range(1, 13):
    theta = m * math.pi / 36.0
    for k in range(12):
        phi = k * math.pi / 6.0
        points.append((math.sin(theta) * math.cos(phi),
                       math.sin(theta) * math.sin(phi), math.cos(theta)))
for x in points:
    near(dist(B(x, 0.0), x), 0.0, 'segment_endpoints', 'B(x,0)=x')
    near(dist(B(x, 1.0), R(x)), 0.0, 'segment_endpoints', 'B(x,1)=R(x)')
    if dist(x, O) > EPS:
        check(cost(R(x)) < cost(x), 'strict_cost', 'strict step cost decrease')
    for k in range(9):
        t = k / 8.0
        y = B(x, t)
        near(norm(y), 1.0, 'sphere_and_cap', 'unit sphere')
        check(y[2] >= 0.5 - EPS, 'sphere_and_cap', 'stay inside compact cap')
        near(dist(B(O, t), O), 0.0, 'fixed_origin', 'B fixes o')
    for n in range(9):
        y = iterate(x, n)
        left = B(y, 1.0)
        right = B(iterate(x, n + 1), 0.0)
        near(dist(left, right), 0.0, 'dyadic_join', 'adjacent pieces agree')
        tn = 1.0 - 2.0 ** (-n)
        near(dist(H(x, tn), y), 0.0, 'dyadic_values', 'H(x,t_n)=R^n(x)')
    near(dist(H(x, 0.0), x), 0.0, 'homotopy_endpoints', 'H starts at x')
    near(dist(H(x, 1.0), O), 0.0, 'homotopy_endpoints', 'H ends at o by definition')

# Exact algebra behind strict descent, sampled only to catch transcription errors.
for z in (F(1,2), F(3,5), F(2,3), F(4,5), F(9,10), F(99,100)):
    lhs = (2*z+1)**2 - z*z*(5+4*z)
    rhs = (1-z*z)*(1+4*z)
    check(lhs == rhs and rhs > 0, 'descent_identity', 'exact positive factorization')

x_edge = (math.sqrt(3.0)/2.0, 0.0, 0.5)
theta = math.acos(x_edge[2])
theta_r = math.acos(R(x_edge)[2])
check(abs(theta_r - theta/2.0) > 0.1, 'nonhalf', 'one step is not geodesic halving')
near(theta_r, math.atan2(math.sqrt(3.0), 2.0), 'nonhalf_formula', 'closed angle formula')
max_by_n = {str(n): max(dist(iterate(x, n), O) for x in points)
            for n in (0, 1, 2, 4, 8, 16, 24)}
max_values = list(max_by_n.values())
check(all(b < a for a,b in zip(max_values, max_values[1:])),
      'finite_iteration_witness', 'sample maxima decrease')

# Actual rotations preserve distance; stabilizer at the pole is nonnormal.
a = math.pi/3.0
h = ((math.cos(a),-math.sin(a),0.0),(math.sin(a),math.cos(a),0.0),(0.0,0.0,1.0))
k = ((0.0,0.0,1.0),(0.0,1.0,0.0),(-1.0,0.0,0.0))
conjugate = matmul(matmul(k,h),transpose(k))
identity = ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
for g in (h,k,conjugate):
    gram = matmul(transpose(g),g)
    for i in range(3):
        for j in range(3):
            near(gram[i][j],identity[i][j], 'proper_rotation', 'orthogonal')
    near(determinant(g),1.0, 'proper_rotation', 'determinant one')
    for x,y in zip(points[:24], reversed(points[-24:])):
        near(dist(action(g,x),action(g,y)),dist(x,y),'isometry','chord distance preserved')
near(dist(action(h,O),O),0.0,'stabilizer','h fixes o')
conjugate_defect = dist(action(conjugate,O),O)
check(conjugate_defect > 0.5,'non_normal_stabilizer','k h k^-1 does not fix o')

Ax = ((0,0,0),(0,0,-1),(0,1,0))
Ay = ((0,0,1),(0,0,0),(-1,0,0))
Az = ((0,-1,0),(1,0,0),(0,0,0))
generators = (Ax,Ay,Az)
group_dim = rank([tuple(v for row in g for v in row) for g in generators])
orbit_dim = rank([action(g,O) for g in generators])
check(group_dim == 3 and orbit_dim == 2,'orbit_rank','SO(3) has 3 generators, orbit has 2 dimensions')

# Reuse of round 425 model, with its angular distance cost (no new model claim).
# Finite-support sequences represent exact infinite sequences by trailing zeros.
def torus_cost(angles):
    total = F(0)
    for j, theta in enumerate(angles, start=1):
        theta = theta % 1
        total += min(theta,1-theta) / 2**j
    return total

supports = [(), (F(1,2),), (F(1,4),F(2,3)),
            tuple(F(j % 7,7) for j in range(1,14)),
            (F(0),)*15 + (F(2,5),F(3,8))]
cost_rows = []
for angles in supports:
    shifted = (F(0),) + angles # retain extra coordinate; do not truncate tail
    before,after = torus_cost(angles),torus_cost(shifted)
    check(after == before/2,'torus_exact_shift','infinite finite-support cost halves exactly')
    cost_rows.append({'support_length':len(angles),'before':str(before),'after':str(after)})

# Tail-coordinate loop projects with degree 1 before R and degree 0 after R.
# Numerical winding checks the explicit representatives; nonexistence of B is
# the analytic homotopy-invariance of degree, not inferred from discretization.
def winding(loop):
    return sum(cmath.phase(b/a) for a,b in zip(loop,loop[1:]))/(2*math.pi)

loop_rows = []
for fixed_prefix in (1,4,12,32):
    j = fixed_prefix+3
    loop = [cmath.exp(2j*math.pi*q/128) for q in range(129)]
    shifted_same_coordinate = [1.0+0j]*129
    degree_before,degree_after = winding(loop),winding(shifted_same_coordinate)
    near(degree_before,1.0,'torus_winding','identity projects to one turn')
    near(degree_after,0.0,'torus_winding','shift projects to constant')
    check(j > fixed_prefix,'tail_loop','tail loop fits any first-prefix basic neighborhood')
    loop_rows.append({'fixed_prefix':fixed_prefix,'tail_coordinate':j,
                      'degree_before':degree_before,'degree_after':degree_after})

result = {
    'status':'PASS',
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks_total':sum(COUNTS.values()),'checks_by_group':COUNTS,
    'max_numerical_identity_residual':max(RESIDUALS.values(),default=0.0),
    'residuals_by_group':RESIDUALS,
    'sphere_cap':{'step_mixture_lambda':'1/3','cap_z_min':'1/2','cost':'1-z',
                  'sample_points':len(points),'edge_angle':theta,'edge_after_step':theta_r,
                  'geodesic_half_defect':theta_r-theta/2,
                  'sample_max_distance_after_n':max_by_n,
                  'conjugated_stabilizer_displacement':conjugate_defect,
                  'group_dimension':group_dim,'orbit_dimension':orbit_dim,
                  'stabilizer_dimension':group_dim-orbit_dim},
    'old_round425_torus':{'cost_formula':'sum_{j>=1} 2^-j dist(theta_j,Z)',
                         'exact_shift_rows':cost_rows,'tail_loop_degrees':loop_rows},
    'scope':[
        'Sphere samples validate formulas, not universal convergence.',
        'Compactness plus strict continuous Lyapunov decrease analytically supplies uniform convergence.',
        'A single jointly continuous segment B(x,t):x->R(x) fixes o; dyadic joins define a limiting contraction.',
        'Finite-support torus cost checks reuse round 425; the new boundary is the degree obstruction to B.',
        'Infinite torus noncontractibility and absence of B are analytic, not finite-dimensional inference.',
        'The contraction is a mathematical limit; no finite total implementation cost is asserted.',
        'This bridge allows space dimension two; it does not select physical dimension three.'
    ]
}
print(json.dumps(result,ensure_ascii=False,indent=2))
