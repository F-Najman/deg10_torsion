from pathlib import Path
import ast,json,re,time,hashlib
from functools import lru_cache
from sage.all import GF,PolynomialRing,PowerSeriesRing,matrix,vector
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C44';B=C/'inputs/inherited_canonical';F=GF(3);start=time.monotonic()
def stream(path):
 text=path.read_text();dec=json.JSONDecoder();i=0
 while i<len(text):
  while i<len(text)and text[i].isspace():i+=1
  if i==len(text):break
  row,i=dec.raw_decode(text,i);yield row
s=(B/'multiplication_22_data.m').read_text();exps=[]
for d,n in[(2,105),(4,245)]:
 v=ast.literal_eval(re.search(r'Std'+str(d)+r'Exponents:=(.*?);',s,re.S)[1]);assert len(v)==n and all(len(e)==36 and sum(e)==d for e in v);exps.append(v)
blocks=[matrix(F,105,245,ast.literal_eval(m[1]))for m in re.finditer(r'Append\(~Tensor22Blocks,Matrix\(GF\(3\),105,245,(\[.*?\])\)\);',s,re.S)];assert len(blocks)==105
for i in range(105):
 for j in range(i):assert blocks[i].row(j)==blocks[j].row(i)
flat=matrix(F,105,105*245,[x for A in blocks for x in A.list()]);monomials=[[(i,e)for i,e in enumerate(row)if e]for row in exps[0]];PF=PolynomialRing(F,'x')
@lru_cache(maxsize=20000)
def eval_jets(serial):
 r=json.loads(serial);d,m=r['degree'],r['multiplicity'];poly=PF(r['field_polynomial']);assert poly.degree()==d and poly.is_monic();K=GF(3**d,name='a',modulus=poly);a=K.gen();pw=[a**i for i in range(d)]
 def dec(v):assert len(v)==d and set(v)<={0,1,2};return sum((K(c)*u for c,u in zip(v,pw)),K(0))
 PS=PowerSeriesRing(K,'t',default_prec=m);zs=[PS([dec(v)for v in row],prec=m)for row in r['canonical_jets']];assert len(zs)==36 and any(z[0]for z in zs)
 powers={}
 for i,z in enumerate(zs):powers[i,1]=z;powers[i,2]=z*z
 rows=[]
 for mon in monomials:
  v=PS.one()
  for i,e in mon:v*=powers[i,e]
  row=[]
  for j in range(m):cc=v[j].polynomial().list();assert len(cc)<=d;row.extend(list(map(int,cc))+[0]*(d-len(cc)))
  rows.append(row)
 return matrix(F,rows)
def replay(name,jets_path,cert_path,expected,survivors,allowed_degrees):
 jj=list(stream(jets_path));assert jj[0]['type']=='header'and jj[-1]['type']=='complete'and jj[-1]['cases']==expected;J={i:{q:[]for q in['T','U']}for i in range(1,expected+1)};dcase={}
 for r in jj:
  if r['type']=='jet':J[r['case']][r['side']].append(r)
  elif r['type']=='case':dcase[r['case']]=r['positive_degree'];assert r['positive_degree']==r['negative_degree']
 if name=='newA5':assert jj[-1]['distinct_source10']==5854 and jj[-1]['distinct_source5']==60
 if name=='residualA7':assert jj[0]['original_cases']==[510,511]
 cert=stream(cert_path);header=next(cert);assert header['type']=='header'and header['cases']==expected;seen=[];ranks=[];done=None
 for r in cert:
  if r['type']=='complete':done=r;break
  assert r['type']=='case';ci=r['case'];assert ci==len(seen)+1;seen.append(ci);T=J[ci]['T'];U=J[ci]['U'];dd=sum(z['degree']*z['multiplicity']for z in T);assert dd in allowed_degrees and sum(z['degree']*z['multiplicity']for z in U)==dd
  if dcase:assert dcase[ci]==dd
  def join(js):
   out=matrix(F,105,0)
   for z in js:out=out.augment(eval_jets(json.dumps({k:z[k]for k in['degree','multiplicity','field_polynomial','canonical_jets']},sort_keys=True)))
   return out
  ET=join(T);EU=join(U);assert ET==matrix(F,105,dd,r['ET2'])and EU==matrix(F,105,dd,r['EU2']);assert ET.rank()==EU.rank()==dd
  f=vector(F,r['F']);h=vector(F,r['H']);assert len(f)==len(h)==105 and f and h and not(f*ET)and not(h*ET)
  MF=matrix(F,105,245,f*flat);MH=matrix(F,105,245,h*flat);assert MF.rank()==105
  Q=matrix(F,245,140,r['quotient_projection']);assert Q.rank()==140 and not(MF*Q)
  A=EU.augment(MH*Q);assert A==matrix(F,105,dd+140,r['matrix']);rr=A.rank();assert rr==r['rank'];assert (rr<105)==(ci in survivors);ranks.append(rr)
  if ci%10==0 or ci in survivors or ci==expected:print('D31A_C44_INDEPENDENT_MATRIX_CASE',name,ci,rr,'seconds',time.monotonic()-start,flush=True)
 assert done and done['survivors']==survivors and len(seen)==expected
 return{'stage':name,'cases':expected,'rank105':sum(r==105 for r in ranks),'survivors':survivors,'degrees':sorted(set(sum(z['degree']*z['multiplicity']for z in J[i]['T'])for i in J)),'jet_sha256':hashlib.sha256(jets_path.read_bytes()).hexdigest(),'matrix_sha256':hashlib.sha256(cert_path.read_bytes()).hexdigest()}
reports=[]
reports.append(replay('newA5',C/'new_finite_jets.normalized.jsonl',C/'new_matrix_results.raw.jsonseq',956,[510,511],{30,60}))
reports.append(replay('residualA7',C/'q7_finite_jets.normalized.jsonl',C/'q7_matrix_results.raw.jsonseq',2,[],{80}))
reports.append(replay('oldA5',B/'sieve_all_inputs.normalized.jsonl',B/'matrix22_all.raw.jsonseq',226,[],{54}))
res={'status':'PASS','stages':reports,'canonical_R2_dimension':105,'canonical_R4_dimension':245,'tensor_full_commutativity_checked':True,'all_F_and_H_vanish_on_positive_restrictions':True,'all_MF_ranks105_and_Q_ranks140':True,'all_selected_product_and_quotient_identities_checked':True,'all_saved_obstruction_matrices_reconstructed_entrywise':True,'cache':str(eval_jets.cache_info()),'elapsed_seconds':time.monotonic()-start}
(C/'independent_matrices_verified.json').write_text(json.dumps(res,indent=2)+'\n');print('D31A_C44_INDEPENDENT_MATRICES_COMPLETE',956,2,226,time.monotonic()-start,flush=True)
