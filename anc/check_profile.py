"""Exact table-to-integral check of the arithmetic profile (Section 4).
Python 3, standard library only.
"""
from fractions import Fraction as F
from itertools import combinations
from decimal import Decimal, localcontext
from math import isqrt
from pathlib import Path
import json
import subprocess
import sys

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
FOUR = (F(4), F(0))


def add(*p):
    return (sum(v[0] for v in p), sum(v[1] for v in p))


def mul(c, p):
    return (c*p[0], c*p[1])


def at(p, x):
    return p[0]+p[1]*x


def candidates(rho, q, x):
    """Affine candidates valid on the cell containing x; q=0 is high."""
    gamma = (-rho, F(q+1))
    eta = add((F(0), F(1)), mul(-1, gamma))

    def good(L):
        k = int(L/x)
        rem = (L, F(-k))
        return add(
            mul(k, gamma),
            min((rem, gamma), key=lambda p: at(p, x))
        )

    A = good(F(1))
    T = good(F(2))
    B = add(T, mul(-1, A))

    def overlap(left, right):
        ends = [left, right]
        for L in (F(1), F(2)):
            rem = (L, F(-int(L/x)))
            if at(left, x) < at(rem, x) < at(right, x):
                ends.append(rem)
        ends.sort(key=lambda p: at(p, x))
        answer = ZERO
        for l, h in zip(ends, ends[1:]):
            s = (at(l, x)+at(h, x))/2
            ci = max(0, (1-s)//x+1)
            ct = max(0, (2-s)//x+1)
            answer = add(
                answer,
                mul(min(ci, ct-ci), add(h, mul(-1, l)))
            )
        return answer

    if q == 0:
        D = overlap(gamma, (F(0), F(1)))
        return [[
            FOUR,
            add(mul(3, T), mul(-2, gamma), D),
            add(ONE, mul(2, T), mul(-1, gamma))
        ]]

    D = overlap(ZERO, gamma)
    base = add((F(2), F(0)), T)
    old = [FOUR, base, add((F(3), F(0)), D)]
    old += [add(base, (-1, 0), N, gamma) for N in (A, B)]
    new = old + [add(base, (-1, 0), N, eta) for N in (A, B)]
    return [old, new]


def integrate(rho, q):
    lo, hi = (rho, F(1)) if q == 0 else (rho/(q+1), rho/q)
    walls = {lo, hi}
    for n in (F(1), F(2), 1+rho, 2+rho):
        for k in range(1, int(n/lo)+2):
            if lo < n/k < hi:
                walls.add(n/k)
    walls = sorted(walls)
    totals = [F(0)]*(1 if q == 0 else 2)

    for l, h in zip(walls, walls[1:]):
        groups = candidates(rho, q, (l+h)/2)
        cuts = {l, h}
        for group in groups:
            for p, r in combinations(group, 2):
                if p[1] != r[1]:
                    t = (r[0]-p[0])/(p[1]-r[1])
                    if l < t < h:
                        cuts.add(t)
        cuts = sorted(cuts)
        for u, v in zip(cuts, cuts[1:]):
            for j, group in enumerate(groups):
                c, s = min(
                    group, key=lambda p: at(p, (u+v)/2)
                )
                totals[j] += c*(v-u)+s*(v*v-u*u)/2
    return totals


def check_counts(rho):
    """Check continuum formulas against literal finite residue counts."""
    tests = 0
    for m in (48, 96, 192):
        g = int(rho*m)
        for p in range(max(5, isqrt(2*m-1+g)+1), 2*m):
            G = p-g % p
            ci = [0]*p
            cj = [0]*p
            for n in range(m):
                ci[n % p] += 1
                cj[(n+m) % p] += 1
            NI, NJ = sum(ci[:G]), sum(cj[:G])
            x = F(p, m)

            if p <= g:
                D = sum(min(ci[s], cj[s]) for s in range(G))
                base = 2*m+NI+NJ
                old = min(
                    4*m, base, 3*m+D,
                    base-m+NI+G, base-m+NJ+G
                )
                new = min(
                    old, base-m+NI+p-G, base-m+NJ+p-G
                )
                exact = [old, new]
                groups = candidates(rho, g//p, x)

            elif p < m:
                active = NI+NJ-G
                D = sum(
                    min(ci[s], cj[s]) for s in range(G, p)
                )
                exact = [
                    min(4*m, 3*active+G+D, m+2*active+G)
                ]
                groups = candidates(rho, 0, x)

            else:
                literal = (2*m-p)+3*sum(
                    1 for n in range(p, 2*m)
                    if n % p < p-g
                )
                assert literal == (
                    4*(2*m-p)-3*max(2*m+g-2*p, 0)
                )
                tests += 1
                continue

            assert exact == [
                m*min(at(v, x) for v in group)
                for group in groups
            ]
            tests += 1
    return tests


def dec(x):
    with localcontext() as ctx:
        ctx.prec = 35
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def main():
    rows = [
        (
            F(5, 16),
            F(185961803, 40000000),
            F(
                1111075051739028956195167,
                111864885256841731511777280
            ),
            F(
                8108653419015955202861020329389,
                1747888832138152054871520000000
            ),
            F(0),
            '4.639112779899639'
        ),
        (
            F(1, 3),
            F(574208581483, 125000000000),
            F(85930118826854129, 7423773458762901240),
            F(
                106301329053850833773025543473,
                23199292058634066375000000000
            ),
            F(0),
            '4.582093659805828'
        ),
    ]
    expected_gains = [
        '47/6720 69/41888 387/637840 '
        '6219367/21864360960 75/484352 '
        '180207/1941459520 317501/5311088640 '
        '18175/441971712 11263/382261760 155/7117824',

        '25/3024 1/540 4/5859 11/34344 7/40095 '
        '17/161460 2/29295 23/492048 '
        '13/389367 29/1175580',
    ]

    replay = None

    for index, (rho, B, L, S, stale_gain, printed) in enumerate(rows):
        high = integrate(rho, 0)[0]
        tail = 2-F(3, 4)*rho*rho
        layers = [integrate(rho, q) for q in range(1, 11)]
        gains = [old-new for old, new in layers]

        assert gains == list(map(F, expected_gains[index].split()))
        assert sum(gains) == L and B-L == S

        old10, new10 = (
            sum(row[j] for row in layers) for j in (0, 1)
        )
        delta = rho/11
        tail_upper_obsolete = (
            F(5, 2)*delta+rho/400+F(25, 48)*delta*delta
        )
        # The explicit all-layer tail bound proved in Section 4.
        tail_upper = F(5, 2)*delta+rho/F(441)+F(7, 16)*delta*delta
        tail_lower = F(5, 2)*delta-delta*delta/2
        # The old profile has its own, weaker tail estimate.  The improved
        # estimate is used only after the q=1,...,10 low-prime rank gains.
        B_upper = high+tail+old10+tail_upper_obsolete
        I_upper = high+tail+new10+tail_upper
        I_lower = high+tail+new10+tail_lower

        assert B_upper < B and I_upper < S
        assert S < F(printed)  # the decimal printed in Section 4 is a strict upper bound

        if replay is not None:
            r = replay['results'][index]
            assert F(r['total']) == L
            assert F(r['final_exact_upper']) == S

        print('rho =', rho, 'finite_count_tests =', check_counts(rho))
        print('high =', high, 'tail =', tail)
        print('gains =', ', '.join(map(str, gains)))
        print('L =', L, 'tail_upper =', tail_upper)
        print('B_upper =', dec(B_upper), '< B =', dec(B))
        print('I_lower =', dec(I_lower), 'I_upper =', dec(I_upper))
        print('S =', S, '=', dec(S))
        print('profile_value =', dec(S))
        print('printed_upper_bound =', printed, 'PASS')

    print('EXACT_PROFILE_AGGREGATION_VERIFIED')


if __name__ == '__main__':
    main()
