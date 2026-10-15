"""Pure-Python interval replay of the two positive endpoint certificates.

This is deliberately independent of python-flint.  It uses mpmath.iv with
exact rational inputs from the frozen density JSON and checks the published
downward endpoint bounds.  It is a cross-check, not a replacement for the
uniform analytic proof.
"""
import json, sys
from fractions import Fraction as Q
from decimal import Decimal
from pathlib import Path
import mpmath as mp
from mpmath import iv

mp.iv.dps = 70

def I(x):
    if isinstance(x, Q):
        return iv.mpf(x.numerator) / x.denominator
    return iv.mpf(x)

def pos_hi(x):
    return x.b

def pos_lo(x):
    return x.a

def chi(x):
    if x == 0:
        return I(0)
    z = I(x)
    return z * (iv.log(abs(z)) - 1)

def psi(x):
    if x == 0:
        return I(0)
    z = I(x)
    return z*z * (iv.log(abs(z))/2 - I(Q(3, 4)))

def rationalize(dens, H):
    den = 1 << 40
    n = [min(max(int(round(Q(v)*den)), 0), 3*den) for v in dens]
    deficit = H*den - sum(n)
    free = [i for i,v in enumerate(n) if 0 < v < 3*den]
    i0 = max(free, key=(lambda i: 3*den-n[i]) if deficit > 0 else (lambda i: n[i]))
    n[i0] += deficit
    assert len(n) == 256 and sum(n) == H*den
    return n, den

def replay(row):
    rho = Q(1, 3); H = 128; N = 256; K = 32
    n, den = rationalize(row['dens'], H)
    h = Q(1, H)
    dens = [Q(v, den) for v in n]
    ell = Q(int(round(Q(row['ell'])*(1 << 30))), 1 << 30)
    lam = iv.log(I(int(round(Q(row['q'])))))
    cells = [(rho+i*h, rho+(i+1)*h, v) for i,v in enumerate(dens) if v]
    L = I(0)
    for i,(a,b,u) in enumerate(cells):
        for j in range(i+1):
            c,d,v = cells[j]
            # The displayed primitive difference is the integral of log|x-y|.
            val = psi(b-c)-psi(a-c)-psi(b-d)+psi(a-d)
            L -= I(u*v*(1 if i == j else 2))*val
    assert pos_lo(L) > 0, ("the signed pair energy L=-int log must be positive", L)
    terms = {Q(-2): Q(1), Q(0): Q(-4), rho: Q(3)}
    for j in range(N+1):
        jump = (dens[j] if j < N else 0) - (dens[j-1] if j else 0)
        c = rho+j*h
        terms[c] = terms.get(c, Q(0))-2*jump
    terms = [(c,v) for c,v in terms.items() if v]
    affine = I(3*rho) * (iv.log(I(rho))-1)
    s = h/K
    penalty = I(0)
    for p in range(N*K):
        a = rho+p*s; b = a+s; mid = (a+b)/2
        la = lam*I(a)+affine; lb = lam*I(b)+affine
        for c,v in terms:
            assert c <= a or c >= b
            if v*(mid-c) > 0:
                f = I(v)*chi(mid-c)
                slope = I(v)*iv.log(abs(I(mid-c)))
                la += f-slope*I(s/2); lb += f+slope*I(s/2)
            else:
                la += I(v)*chi(a-c); lb += I(v)*chi(b-c)
        ua = pos_hi(I(ell)-la); ub = pos_hi(I(ell)-lb)
        if ua <= 0 and ub <= 0:
            continue
        if ua >= 0 and ub >= 0:
            v = I(s)*(ua+ub)/2
        else:
            pos,neg = (ua,ub) if ua > 0 else (ub,ua)
            v = I(s)*pos*pos/(2*(pos-neg))
        penalty += pos_hi(v)
    T = rho+N*h
    V = lambda t: chi(t+2)-4*chi(t)+3*chi(t-rho)+lam*I(t)+affine
    assert pos_lo(lam*I(T-rho)) > I(3*rho+2)
    assert pos_lo(V(T)-2*iv.log(abs(I(T-rho)))) > I(ell)
    prim = lambda t: psi(t+2)-4*psi(t)+3*psi(t-rho)+lam*I(t*t)/2+affine*I(t)
    integralV = I(0)
    for a,b,v in cells:
        integralV += I(v)*(prim(b)-prim(a))
    upper = I(Q(3,2))+integralV+L
    lower = I(Q(3,2))-L+I(ell)-3*penalty
    target = I(Q('4.9972574688' if int(round(row['q'])) == 27 else '5.0224702721'))
    assert pos_hi(lower) <= pos_lo(upper), ("dual lower bound exceeds E(sigma)", lower, upper)
    assert pos_lo(lower) >= target, ("directed lower bound misses printed target", lower, target)
    print(json.dumps({'q': int(round(Q(row['q']))), 'lower': str(lower),
                      'lower_endpoint': str(pos_lo(lower)), 'upper': str(upper),
                      'pair_energy': str(L), 'target': str(target), 'pass': True}))
    return {'lower': str(lower), 'upper': str(upper),
            'pair_energy': str(L), 'target': str(target), 'pass': True}

if __name__ == '__main__':
    data = json.loads(Path(sys.argv[1]).read_text(), parse_float=Decimal)
    for row in data:
        replay(row)
