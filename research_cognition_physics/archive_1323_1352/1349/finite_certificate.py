"""Exact rational enclosure of four nonseparability witnesses, including two per physical group."""
from fractions import Fraction as F
from math import factorial


def y_amplitude(shift, time):
    # Real Taylor coefficients for H^k|x>, phases recorded by k modulo four.
    x,y = F(1),F(0)
    re,im = F(0),F(0)
    for k in range(1,81):
        x,y = y,x+shift*y
        x,y = x*time/k,y*time/k
        if k%2 == 0:
            re += y if k%4==0 else -y
        else:
            im += y if k%4==3 else -y
    return re,im


def prob(z, sign, time):
    delta = F(13,10)*(z-1)
    mass = F(7,10)*(2+sign)
    a,b = y_amplitude(delta+mass,time),y_amplitude(delta-mass,time)
    return sum(((u-v)/2)**2 for u,v in zip(a,b))


def records(channels,time):
    symmetric = sum(w*prob(z,1,time) for z,w in channels)
    antisymmetric = sum(w*prob(z,-1,time) for z,w in channels)
    return symmetric,(symmetric+antisymmetric)/2


# Operator infinity bound (row sum) on the Hermitian 2x2 generators.
error = F(13)**81/factorial(81)  # T*||H|| < 13 in all cases
assert error < F(1,10**27)
prob_error = 2*error+error**2
det_error = 4*prob_error+2*prob_error**2
groups = [[(F(1),F(1,3)),(F(5,2),F(2,3))],
          [(F(0),F(1,9)),(F(1),F(1,3)),(F(3),F(5,9))]]
times = [F(7,5),F(11,5)]
table = [[records(g,t) for t in times] for g in groups]
det = lambda a,b: a[0]*b[1]-a[1]*b[0]
for group in table:
    a,b=group
    value=det(a,b)
    assert value-det_error > F(1,200)
    print("time and bank determinant",float(value),
          "held-out minus product estimate",float(b[1]-b[0]*a[1]/a[0]))
for i in range(2):
    value=det(table[0][i],table[1][i])
    assert value-det_error > F(1,10000)
    print("group and bank determinant",float(value))
