from pathlib import Path
from functools import lru_cache
import ast,re,json,time,hashlib
from sage.all import GF,PolynomialRing,PowerSeriesRing,matrix
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C44';start=time.monotonic();k=GF(3);PF=PolynomialRing(k,'y')
s=(C/'inputs/canonical_data.m').read_text();assert 'AFBaseCoordinate:=2;'in s
get=lambda n:ast.literal_eval(re.search(n+r':=(.*?);',s,re.S)[1].replace('<','[').replace('>',']'))
af=get('AFPolynomialData');canon=get('CanonicalCoordinateData');assert len(af)==20 and len(canon)==36 and all(len(row)==19 for row in canon);assert all(p[1]==[1]for p in af)and all(p[1]==[1]for row in canon for p in row)
polys=[PF(q[0])for q in af];derivs=[q.derivative()for q in polys]
mons=sorted({(j,l)for row in canon for j,p in enumerate(row)for l,a in enumerate(p[0])if a});A=matrix(k,[[row[j][0][l]if l<len(row[j][0])else 0 for j,l in mons]for row in canon]);assert A.nrows()==36
count=0
@lru_cache(maxsize=20000)
def verify(serial):
 global count
 r=json.loads(serial);d,m=r['degree'],r['multiplicity'];poly=PF(r['field_polynomial']);assert poly.degree()==d;L=GF(3**d,'a',modulus=poly);a=L.gen()
 def dec(v):assert len(v)==d and set(v)<={0,1,2};return PF(v)(a)
 PS=PowerSeriesRing(L,'T',default_prec=m);T=PS.gen();xy0=[dec(v)for v in r['xy']];xy=[PS([dec(v)for v in row],prec=m)for row in r['xy_jets']];assert len(xy)==2 and all(v[0]==w for v,w in zip(xy,xy0));x,y=xy;x0,y0=xy0
 xp=[PS.one()];yp=[PS.one()]
 for j in range(1,20):xp.append(xp[-1]*x)
 for j in range(1,max(l for j,l in mons)+1):yp.append(yp[-1]*y)
 fx=sum((j*polys[j](y0)*x0**(j-1)for j in range(1,len(polys))),L(0));fy=sum((derivs[j](y0)*x0**j for j in range(len(polys))),L(0));dir=r['uniformizer_coordinate']
 assert (dir==1 and fy!=0)or(dir==2 and fy==0 and fx!=0)
 assert (x==PS(x0)+T if dir==1 else y==PS(y0)+T)
 ff=sum((polys[j](y)*xp[j]for j in range(len(polys))),PS.zero());assert ff.is_zero()
 rows=[]
 for j,l in mons:
  v=xp[j]*yp[l];row=[]
  for h in range(m):cc=v[h].polynomial().list();assert len(cc)<=d;row.extend(list(map(int,cc))+[0]*(d-len(cc)))
  rows.append(row)
 MM=A*matrix(k,rows);raw=[]
 for row in MM.rows():raw.append(PS([sum((L(row[d*h+i])*a**i for i in range(d)),L(0))for h in range(m)],prec=m))
 nz=next(j for j,v in enumerate(raw)if v[0]);unit=raw[nz];zs=[v/unit for v in raw];expected=[PS([dec(v)for v in row],prec=m)for row in r['canonical_jets']];assert zs==expected and len(zs)==36
 count+=1
 if count%25==0:print('D31A_C44_JET_GEOMETRY_UNIQUE',count,'degree',d,'precision',m,'seconds',time.monotonic()-start,flush=True)
 return True
reports=[]
for name,path in [('newA5',C/'new_finite_jets.normalized.jsonl'),('residualA7',C/'q7_finite_jets.normalized.jsonl'),('oldA5',C/'inputs/inherited_canonical/sieve_all_inputs.normalized.jsonl')]:
 nr=0
 for line in path.read_text().splitlines():
  if not line.strip():continue
  r=json.loads(line)
  if r['type']!='jet':continue
  verify(json.dumps({key:r[key]for key in['degree','multiplicity','field_polynomial','xy','xy_jets','canonical_jets','uniformizer_coordinate']},sort_keys=True));nr+=1
 reports.append({'stage':name,'all_jet_blocks_verified':nr,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()});print('D31A_C44_JET_GEOMETRY_STAGE_COMPLETE',name,nr,time.monotonic()-start,flush=True)
res={'status':'PASS','stages':reports,'unique_jet_geometries':count,'all_affine_defining_equations_verified_to_full_saved_precision':True,'all_uniformizer_parameters_and_nonsingular_partials_verified':True,'all36_global_canonical_coordinate_functions_verified_at_every_saved_coefficient':True,'normalization':'first nonzero constant canonical coordinate','coordinate_model':'root x, base y, AF degree19, genus36','cache':str(verify.cache_info()),'elapsed_seconds':time.monotonic()-start}
(C/'independent_jet_geometry_verified.json').write_text(json.dumps(res,indent=2)+'\n');print('D31A_C44_INDEPENDENT_JET_GEOMETRY_COMPLETE',count,time.monotonic()-start,flush=True)
