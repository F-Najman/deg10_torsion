from pathlib import Path
import ast,re,json,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,matrix,vector
H=Path(__file__).resolve().parent;k=GF(3);PR=PolynomialRing(k,'a');tt=time.monotonic()
def lit(name):
 s=(H/'quotient_replay_points.m').read_text();v=re.search(r'\b'+name+r':=(.*?);',s,re.S).group(1)
 v=re.sub(r'\[\s*\*','[',v);v=re.sub(r'\*\s*\]',']',v);v=v.replace('<','(').replace('>',')').replace('true','True').replace('false','False');return ast.literal_eval(v)
pts=lit('SourceTargetPoints');cs=lit('CaseClasses');assert len(pts)==676 and len(cs)==1023
R2=PolynomialRing(k,names=('x','y'));vx,vy=R2.gens();polys={}
for n in[23,46]:polys[n]=eval('('+(H/'inputs'/('FFFc'+str(n)+'.txt')).read_text().replace('^','**')+')',{'__builtins__':{}},{'x':vx,'y':vy})
def xy(E,P):
 a1,a2,a3,a4,a6=E.ainvs();x,y=P.xy();v3=2*y+a3+a1*x;v2=3*x+a2;v4=3*x*x+2*x*a2+a4-a1*y
 z1=2*v4/v3+a1;z2=v2-(v4/v3)**2-a1*v4/v3;b=-z2**3/v3**2;c=-(z1*z2-v3)/v3;r=b/c;s=c*c/(b-c);v=r*s-2*r+1
 return (s-r)/v,v/(s*s-s-r+1)
def key(x,y,d):
 for h in[x,y,x+y,x-y]:
  if h.minpoly().degree()==d:break
 else:raise AssertionError('primitive')
 def co(v):a=v.polynomial().list();return list(a)+[0]*(d-len(a))
 A=matrix(k,[co(h**i)for i in range(d)]);AI=A.inverse()
 return(tuple(map(int,h.minpoly().list())),tuple(map(int,vector(k,co(x))*AI)),tuple(map(int,vector(k,co(y))*AI)))
fields={};units={u for u in range(1,23)if u%2==1};source10=set();sourceother={d:set()for d in[4,5,6]};candidates=set(range(1,23));fullclasses={}
for i,a in enumerate(pts,1):
 d,cp,co,ai,dr,norm,txy,cl,qv=a;cp=tuple(cp)
 if cp not in fields:assert PR(cp).is_irreducible();fields[cp]=GF(3**d,'a',modulus=PR(cp))
 L=fields[cp];u=L.gen();vv=[PR(z)(u)for z in co];assert vv[2]!=0;x,y=vv[0]/vv[2],vv[1]/vv[2];assert polys[46](x,y)==0
 r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1);E=EllipticCurve(L,[1-c,-b,-b,0,0]);P=E(0,0)
 assert list(E.ainvs())==[PR(v)(u)for v in ai]and 46*P==E(0)and 23*P!=E(0)and 2*P!=E(0)
 xx,yy=xy(E,2*P);assert (xx,yy)==tuple(PR(v)(u)for v in txy)and polys[23](xx,yy)==0
 actual=min(h for h in range(1,d+1)if d%h==0 and xx**(3**h)==xx and yy**(3**h)==yy);assert actual==dr and norm*dr==d
 assert len(cl)==len(qv)==2 and cl[1]==qv[1]==0
 candidates={v for v in candidates if (v*norm*cl[0]-qv[0])%23==0};assert candidates
 # Reconstruct all marked source diamond keys; fields may have arbitrary primitive polynomials.
 orbit={key(*xy(E,a*P),d)for a in units};assert len(orbit)==11
 if d==10:assert not(source10&orbit);source10|=orbit
 else:sourceother[d]|=orbit
 if i%100==0:print('D31A_C46_INDEPENDENT_BINDING',i,round(time.monotonic()-tt,2),flush=True)
assert len(source10)==6292 and {d:len(v)for d,v in sourceother.items()}=={4:55,5:44,6:55}
assert len(candidates)==1
original=json.loads((H/'source_cases.json').read_text())['cases'];lookup={(tuple(a['support_type']),a['index']):a for a in original};assert len(lookup)==1023
exc=[]
for typ,idx,refs,cl,reject in cs:
 value=sum(mult*pts[i-1][8][0]for i,mult in refs)%23;assert cl==[value,0]and reject==(value!=0)
 old=lookup[tuple(typ),idx];assert old['quotient_reject']==reject
 if not reject:exc.append([typ,idx])
assert len(exc)==37
out=dict(status='PASS',source_components=676,source_degree10_points=6292,source_small_counts={d:len(v)for d,v in sourceother.items()},support_cases=1023,quotient_rejections=986,exceptional_cases=exc,
 quotient_generator_orientation=list(candidates)[0],all_marked_source_and_target_points_verified=True,all_residue_norms_and_class_sums_verified=True,seconds=time.monotonic()-tt)
(H/'independent_quotient_arithmetic.json').write_text(json.dumps(out,indent=2)+chr(10));print('D31A_C46_INDEPENDENT_QUOTIENT_COMPLETE',out,flush=True)
