from pathlib import Path
import ast,re,json,time,math
from sage.all import GF,PolynomialRing,EllipticCurve
P=Path(__file__).resolve().parent;t0=time.monotonic()
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']').replace('true','True').replace('false','False').strip())
def assignment(path,name):
 s=(P/path).read_text();return lit(re.search(name+r':=(.*?);',s,re.S)[1])
pol=assignment('matched_levels.m','FieldPolynomial');R=PolynomialRing(GF(2),'a');F=GF(2**10,'a',modulus=R(pol));a=F.gen()
def dec(c):return sum(F(v)*a**i for i,v in enumerate(c))
def code(v):return list(map(int,v._vector_()))
units=[u for u in range(1,57)if math.gcd(u,57)==1];assert len(units)==36
assert {sgn*pow(5,k,57)%57 for sgn in [1,-1]for k in range(18)}==set(units)
def key(Q):return min(code((u*Q)[0])+code((u*Q)[1])for u in units)
raw=assignment('matched_levels.m','AllSourceMatches');matched=assignment('matched_levels.m','MatchedLevels');missing=assignment('matched_levels.m','MissingLevels');previous=assignment('selected_levels.m','SelectedLevels')
assert len(raw)==12 and len(matched)==6 and len(missing)==2 and len(previous)==6
filterdata=lit((P/'noncuspidal_filter_v2_data.txt').read_text());residual={tuple(v[3])for v in filterdata[1]};assert len(residual)==3
byj={}
for item in raw:
 rec,jdual,needed,old=item;idx,jp,iv,aa,pp,k=rec;E=EllipticCurve(F,[dec(v)for v in aa]);Q=E(*[dec(v)for v in pp]);assert Q.order()==57 and list(E.j_invariant().minpoly())==jp and key(Q)==k
 byj.setdefault(tuple(jp),[]).append((rec,E,Q))
assert len(byj)==6
orbit_counts=[]
for jp,rows in byj.items():
 rec,E,Q=rows[0];points=E.points();assert len(points)==1026
 pts=[Q for Q in points if Q and 57*Q==E(0)and 19*Q!=E(0)and 3*Q!=E(0)]
 keys={tuple(key(Q))for Q in pts};native={tuple(row[0][5])for row in rows};assert keys==native and len(pts)==36*len(keys)
 orbit_counts.append(dict(jpoly=list(jp),primitive_points=len(pts),cyclic_subgroup_orbits=len(keys)))
 print('C57_DECIC_ALL_PHYSICAL_CYCLIC_SUBGROUPS_PASS',jp,len(keys),flush=True)
required=[];matches=[]
for rec,jdual,needed,old in raw:
 idx,jp,iv,aa,pp,k=rec;E=EllipticCurve(F,[dec(v)for v in aa]);Q=E(*[dec(v)for v in pp]);f3=E.isogeny(19*Q);assert f3.degree()==3 and f3(Q).order()==19;f19=f3.codomain().isogeny(f3(Q));assert f19.degree()==19 and f19(f3(Q))==f19.codomain()(0)
 dpoly=list(map(int,f19.codomain().j_invariant().minpoly()));assert dpoly==jdual and needed==(tuple(dpoly)in residual)
 testids=[i+1 for i,v in enumerate(previous+missing)if v[1]==jp and v[5]==k];assert len(testids)<=1
 if needed:assert len(testids)==1;required+=testids
 matches.append(dict(source_jpoly=jp,dual_jpoly=dpoly,passes_containing_j_filter=needed,test_ids=testids))
 print('C57_DECIC_INDEPENDENT_3_THEN_19_DUAL_MATCH_PASS',idx,needed,testids,flush=True)
assert len(required)==6 and len(set(required))==6
out=dict(complete=True,field_degree=10,all_physical_j_orbits=6,all_physical_subgroup_orbits=12,residual_dual_j_orbits=3,required_physical_subgroup_orbits=6,required_test_indices=sorted(required),retained_tests=8,all_points_enumerated=True,independent_dual_j_by_3_then_19_isogenies=True,orbit_counts=orbit_counts,matches=matches)
(P/'source_coverage_verified.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_DECIC_SOURCE_COVERAGE_INDEPENDENT_SAGE_PASS',sorted(required),'seconds',time.monotonic()-t0,flush=True)
