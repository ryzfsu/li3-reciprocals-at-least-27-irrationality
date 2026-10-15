"""Pure-Python interval replay of the two negative signed-contour certificates.

Uses only mpmath.iv, Fraction, and the frozen rational masses/boxes.  It
reconstructs pair energies and checks every archived interval box against its
downward-rounded target; it does not import python-flint or SciPy.
"""
import json, sys
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
from mpmath import iv

iv.dps = 60
RHO = Q(5, 16)

def I(x):
    if isinstance(x, Q):
        return iv.mpf(x.numerator) / x.denominator
    return iv.mpf(x)

def CI(re, im=0):
    return iv.mpc(I(re), I(im))

def neg_hi(x):
    return abs(x).b

def neg_lo(x):
    return x.a

def zlog(z):
    return iv.log(z)

def zg(z):
    if z == 0:
        return CI(0)
    return z*zlog(z)-z

def zf(z):
    if z == 0:
        return CI(0)
    return z*z*(zlog(z)/2-I(Q(3,4)))

def rf(t):
    if t == 0:
        return I(0)
    x=I(t)
    return x*x*(iv.log(abs(x))/2-I(Q(3,4)))

def rj(t):
    if t == 0:
        return I(0)
    x=I(t)
    return x*(iv.log(abs(x))-1)

def is_positive(x):
    return neg_lo(x) > 0

def replay(q, root):
    data=json.loads((root/'certificate_data.json').read_text())
    record=data['negative'][str(q)]
    counts=[int(x) for x in record['counts']]
    den=int(record['denominator'])
    assert len(counts)==256 and sum(counts)==den and all(x>=0 for x in counts)
    edges=[Q(3*j**3,256**3) for j in range(257)]
    cells=[]; D={}
    for j,k in enumerate(counts):
        if k:
            a,b=edges[j],edges[j+1]
            p=Q(k,den); slope=p/(b-a)
            cells.append((a,b,p,slope))
            D[a]=D.get(a,Q(0))+slope
            D[b]=D.get(b,Q(0))-slope
    D={e:v for e,v in D.items() if v}
    nodes=sorted(D)
    Ac=iv.mpc(I(4)/iv.sqrt(I(17)),I(1)/iv.sqrt(I(17)))
    beta=iv.mpc(-I(4)/iv.sqrt(I(17)),I(1)/iv.sqrt(I(17)))
    Cc=iv.mpc(I(Q(-15,17)),I(Q(-8,17)))
    lam=iv.log(I(q)); linear=lam*Ac.real+iv.pi*Ac.imag

    def zval(x):
        return I(RHO) + Ac*I(x)

    def potential(u):
        z=zval(u)
        zz=lambda w: zg(w)+w
        return (zz(z+2)-zz(z)-2).real \
            -3*((zz(z)-zz(z-I(RHO))).real-I(RHO)*iv.log(I(RHO))) \
            +lam*z.real+iv.pi*z.imag

    def potential_primitive(uq):
        u=I(uq); z=zval(u)
        b_primitive=(zf(z+2)-zf(z))/Ac
        h_primitive=(zf(z)-zf(z-I(RHO)))/Ac \
            +I(RHO)*(1-iv.log(I(RHO)))*u
        return (b_primitive-3*h_primitive+lam*I(RHO)*u
                +linear*u*u/2).real

    def Wpoint(uq):
        u=I(uq); out=potential(u)
        for a,b,p,slope in cells:
            aa,bb=I(a),I(b)
            same=rj(uq-a)-rj(uq-b)
            opp=((zg(Ac*u+beta*bb)-zg(Ac*u+beta*aa))/beta).real
            out-=I(slope)*(same+opp)
        return out

    # Exact cell primitives give the pair-log energy, with interval logs.
    L_sigma=I(0)
    for ii,(a,b,p,slope) in enumerate(cells):
        aa,bb=I(a),I(b)
        for jj in range(ii+1):
            c,e,pp,other=cells[jj]
            cc,ee=I(c),I(e)
            # This corner difference is already -int int log|u-v|.
            same=rf(b-e)-rf(a-e)-rf(b-c)+rf(a-c)
            # The signed energy uses -log|z-w| on both rays.
            opp=-((zf(Ac*bb+beta*ee)-zf(Ac*aa+beta*ee)
                   -zf(Ac*bb+beta*cc)+zf(Ac*aa+beta*cc))/(Ac*beta)).real
            term=I(slope*other)*(same+opp)/2
            L_sigma += term if ii==jj else 2*term
    assert neg_lo(L_sigma)>0, ("the signed pair energy L=-int log must be positive", L_sigma)

    integralV=I(0)
    for a,b,p,slope in cells:
        integralV += I(p/(b-a))*(potential_primitive(b)-potential_primitive(a))
    upper=I(Q(3,2))+integralV+L_sigma

    def derivative_residual(u, omit=None):
        z=I(RHO)+Ac*u
        out=(Ac*(zlog(z+2)-4*zlog(z)+3*zlog(Ac))).real+linear
        cusp=(omit==Q(0))
        if not cusp:
            out += 3*Ac.real*iv.log(u)
        for e,slope_q in D.items():
            slope=I(slope_q)
            if e!=omit:
                out -= slope*iv.log(abs(u-I(e)))
            if cusp and e==0:
                out += slope*(Cc*zlog(Ac)).real
            else:
                out += slope*(Cc*zlog(Ac*u+beta*I(e))).real
        return out

    boxes=record['boxes']
    boxes=sorted(boxes,key=lambda b:Q(b['l']))
    cursor=Q(0)
    target=I('4.972' if q==27 else '4.990')
    passed=0; minimum=None; minimum_iv=None
    for box in boxes:
        l,r=Q(box['l']),Q(box['r'])
        node=None if box['node'] is None else Q(box['node'])
        assert l==cursor and r>l; cursor=r
        has=node is not None
        if has:
            assert node in D and l<=node<=r
        for e in D:
            if l<=e<=r:
                assert has and e==node
        # Form the hull from directed interval enclosures of exact fractions.
        ul=I(l); ur=I(r)
        u=iv.mpf([ul.a,ur.b])
        if has:
            delta=max(node-l,r-node); de=I(delta)
            residual=derivative_residual(u,node)
            if node==0:
                coeff=3*Ac.real-I(D.get(Q(0),0))+I(D.get(Q(0),0))*Cc.real
            else:
                coeff=-I(D[node])
            bound=Wpoint(node)-de*neg_hi(residual)-neg_hi(coeff)*de*(1+neg_hi(iv.log(de)))
        else:
            mid=(l+r)/2; de=I((r-l)/2)
            residual=derivative_residual(u,None)
            bound=Wpoint(mid)-de*neg_hi(residual)
        full=I(Q(3,2))-L_sigma+bound
        assert neg_lo(full) > target, (q, l, r, node, full)
        fl=neg_lo(full)
        if minimum is None or fl < minimum:
            minimum=fl
            minimum_iv=full
        passed+=1
    assert cursor==3 and passed==len(boxes)
    support=cells[-1][1]
    tail=linear-I(RHO)-2/(3-I(support))
    w3=I(Q(3,2))-L_sigma+Wpoint(Q(3))
    assert is_positive(tail) and neg_lo(w3)>target
    assert minimum_iv is not None
    assert neg_hi(minimum_iv)<=neg_lo(upper), ("dual lower bound exceeds E(sigma)", minimum_iv, upper)
    assert neg_lo(minimum_iv)>=target, ("directed lower bound misses printed target", minimum_iv, target)
    receipt={'lower':str(minimum_iv),'upper':str(upper),'pair_energy':str(L_sigma),
             'minimum_box_lower':str(minimum),'tail_lower':str(neg_lo(tail)),
             'endpoint_value_lower':str(neg_lo(w3)),'target':str(target),
             'boxes':passed,'pass':True}
    print(json.dumps({'q':q,'boxes':passed,'minimum_box_lower':str(minimum),
                      'pair_energy':str(L_sigma),'lower':str(minimum_iv),
                      'upper':str(upper),'tail_lower':str(neg_lo(tail)),
                      'endpoint_value_lower':str(neg_lo(w3)),'target':str(target),
                      'pass':True}),flush=True)
    return receipt

if __name__=='__main__':
    root=Path(sys.argv[1])
    replay(27,root)
    replay(28,root)
