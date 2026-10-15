"""Exact determinant identities by degree-bounded interpolation; no external files."""
import json, math, time, platform, re
from pathlib import Path
import flint
from flint import fmpq_poly as QP, fmpz_poly as ZP, fmpz_mat
ROOT=Path(__file__).resolve().parent
x=QP([0,1])
receipts=[]
for m,g,name,e,f,prime in [(12,4,'H211',12,83,29),(16,5,'H378',24,150,37)]:
 start=time.time();B=3*g-3;w=QP([1])
 for s in range(1,g+1):w*=x-s;w/=s
 w=w**3;M=[];scales=[]
 for i in range(m):
  row=[]
  for j in range(m,2*m):
   P,rem=divmod(w,(x+i)*(x+j))
   ci,cj=w(-i),w(-j);di,dj=ci/(j-i),-cj/(j-i)
   assert rem==di*(x+j)+dj*(x+i)
   der=P.derivative().derivative(); assert der.degree()<=B-1
   vals=[der(s) for s in range(1,B+1)];ds=[]
   while vals:
    ds.append(vals[0]); vals=[vals[t+1]-vals[t] for t in range(len(vals)-1)]
   ent=sum((c*(x-1)**(B-1-t) for t,c in enumerate(ds)),QP([]))
   for n,d in [(i,di),(j,dj)]:
    ent-=2*d*sum((x**(n-s)/s**3 for s in range(1,n+1)),QP([]))*(x-1)**B
   assert ent.degree()<=j+B-1
   row.append(ent)
  scale=math.lcm(*(int(c.q) for v in row for c in v.coeffs()))
  scales.append(scale);M.append([ZP([int(c*scale) for c in v.coeffs()]) for v in row])
 rec=json.loads((ROOT/(name+'_exact_seed.json')).read_text());H=ZP(rec['coefficients_ascending'])
 D=int(rec['C_denominator']);S=math.prod(scales)
 assert rec['C_numerator']==1
 assert D==math.prod(int(p)**z for p,z in rec['D_factorization'].items())
 assert math.gcd(*rec['coefficients_ascending'])==1
 assert rec['a_factor']==e
 assert rec.get('am1_numerator_factor',rec.get('am1_factor'))==f
 assert rec['row_scales']==scales
 deg_bound=sum(j+B-1 for j in range(m,2*m))
 assert H.degree()+e+f==deg_bound
 # Both cleared sides have degree <= deg_bound. Exact agreement at
 # deg_bound+1 distinct integers proves equality in Q[x].
 lo=-(deg_bound//2)
 for count,t in enumerate(range(lo,lo+deg_bound+1),1):
  det=int(fmpz_mat([[int(v(t)) for v in row] for row in M]).det())
  assert D*det==S*t**e*(t-1)**f*int(H(t)),(name,t)
  if count%100==0:print(name,'exact points',count,flush=True)
 vals=[int(H(t))%prime for t in range(prime)]
 assert all(vals)
 assert vals==rec.get('no_root_values',rec.get('mod37_values'))
 out=dict(seed=[m,g],degree_bound=deg_bound,exact_points=deg_bound+1,interval=[lo,lo+deg_bound],a_factor=e,am1_den_power=m*B-f,root_exclusion_prime=prime,root_values=vals)
 receipts.append(out);print(json.dumps(out),flush=True)
tex = ''.join((ROOT.parent / name).read_text()
              for name in ('main.tex', 'body.tex', 'appendices.tex'))
for degree, prime, key in ((211, 29, 'no_root_values'),
                           (378, 37, 'mod37_values')):
    # The appendix uses a gathered display, while the body mentions H_degree
    # several times before the actual residue list.  Anchor at the seed label
    # and capture only the parenthesized vector inside that display.
    seed_label = 'For the positive seed' if degree == 211 else 'For the negative seed'
    match = re.search(
        rf'{re.escape(seed_label)}.*?values modulo \$?{prime}\$?.*?'
        rf'\\begin\{{gathered\}}(.*?)\\end\{{gathered\}}',
        tex, re.S)
    assert match, f'missing displayed residues for H_{degree}'
    shown = [int(x) for x in re.findall(r'-?\d+', match.group(1))]
    rec = next(r for r in receipts if r['root_exclusion_prime'] == prime)
    assert shown == rec['root_values']

(ROOT/'seed_identity_receipt.json').write_text(
    json.dumps(dict(python=platform.python_version(),
                    python_flint=flint.__version__,
                    identities=receipts), indent=2)+'\n')
print('seed identities=PASS exact_points=307/553 residues=PASS')
