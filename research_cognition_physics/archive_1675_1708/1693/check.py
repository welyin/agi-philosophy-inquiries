"""Inventory the two flavor actions on the same fermions; test minimal repairs and graph constraints."""
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def signed(q):
    return ((q+1) % 3)-1


def cyclic_inside(charges, n):
    if n == 3:
        return sum(signed(q)**3 for q in charges) % 9 == 0
    if n == 12:
        return sum(q**3 for q in charges) % 36 == 0 and sum(charges) % 6 == 0
    raise ValueError(n)


def vertices(cfg):
    good = []
    for values in itertools.product(range(3),repeat=8):
        ten, five = values[:3], values[3:6]
        up, down = values[6:]
        if (up+down) % 3 or any((2*t+up) % 3 for t in ten):
            continue
        if any((ten[i]+five[j]+down) % 3 for i,j in cfg['down_edges']):
            continue
        for n in (1, 2):
            if all((l+up+n) % 3 == 0 for l in five):
                good.append([*ten,*five,up,down,n])
    assert len(good) == 2
    return good


def run():
    inp = read('inputs.json')
    sec = inp['auxiliary_3_signed']
    old, flav = read(inp['ordinary_12']),read(inp['old_mediators'])['flavor_charge']
    assert inp['neutral_charges_12']['N'] == old['neutral_pair_example']
    fields = []
    for q in old['su5_generation_charges_10']:
        fields += [(q,sec['ten'])]*10
    for q in old['su5_generation_charges_5bar']:
        fields += [(q,sec['fivebar'])]*5
    for field in ['P','Pbar','T','Tbar']:
        fields += [(flav[field],0)]*5
    fields += [(flav[f],0) for f in old['neutral_charged_under_flavor']]
    # Higgsinos are paired conjugates. Other neutral multipliers and real gauginos
    # are trivial for the finite factors, hence can be omitted from charge sums.
    fields += [(0,sec['up_higgs']),(0,sec['down_higgs'])]*2
    for k in ['N','H']:
        fields += [(q,sec[k]) for q in inp['neutral_charges_12'][k]]
    fields += [(inp['neutral_charges_12']['B'],sec['B'])]
    aux = [x for _,x in fields]
    assert sum(aux) == sum(x**3 for x in aux) == -31
    solutions = []
    for count in range(5):
        repairs = [list(row) for row in itertools.combinations_with_replacement((-1,0,1),count)
                   if cyclic_inside(aux+list(row),3)]
        solutions.append(repairs)
    assert solutions == [[],[],[],[],[[1,1,1,1]]]
    repaired = fields+[(0,1)]*4
    mixed = Counter((signed(q),x) for q,x in repaired if signed(q) and x)
    assert mixed == {(1,-1):10,(-1,1):1}
    # Check every ordinary 3- and 12-element cyclic restriction separately.
    # Completeness for non-cyclic bundles is proved by the character splitting.
    for power3, poweraux in itertools.product(range(3), repeat=2):
        assert cyclic_inside([power3*q+poweraux*x for q,x in repaired],3)
    for power12, poweraux in itertools.product(range(12),range(3)):
        assert cyclic_inside([power12*q+4*poweraux*x for q,x in repaired],12)
    index_aux = 3*sum([sec['ten']]*3)+sum([sec['fivebar']]*3)
    assert index_aux == -9 and index_aux % 3 == 0
    # The extra odd singlets and N0 share all bilinear selection charges.
    b = sec['B']
    degrees = [k for k in range(5) if (2+k*b)%3 == 0]
    assert degrees == [2]
    return dict(allowed_family_higgs_n_assignments=vertices(inp),
                auxiliary_signed_linear_and_cubic=[sum(aux),sum(x**3 for x in aux)],
                repairs_0_to_4=solutions, mixed_3primary_multiplicities={
                    ','.join(map(str,k)):v for k,v in mixed.items()},
                cubic_after_repair=sum(x**3 for _,x in repaired),
                mixed_SU5_index_new=index_aux,
                tested_cyclic_restrictions={'Z3_maps':9,'Z12_maps':36},
                first_B_degree_in_spectator_mass_and_competing_N0_contact=degrees[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    answer = run()
    out = HERE/'results.json'
    if parser.parse_args().check:
        assert answer == json.loads(out.read_bytes())
        print('PASS: connected weak texture, four-singlet minimum, finite product and mixed characters')
    else:
        with out.open('x',encoding='utf8') as fp:
            fp.write(json.dumps(answer,indent=2)+'\n')
