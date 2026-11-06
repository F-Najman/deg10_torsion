from pathlib import Path
import json,hashlib,sys
P=Path(__file__).parent
qs=list(map(int,sys.argv[1:]))or[7,13]
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
def mag(v):
 if isinstance(v,int):return str(v)
 f=all(isinstance(z,int)for z in v);return('['if f else'[*')+','.join(mag(z)for z in v)+(']'if f else'*]')
a=rd('low_hecke_orbits.m');out=[]
for q,n,pure in[(7,2,True),(13,4,False)]:
 if q not in qs:continue
 points={};cases={};models=[]
 for s in range(n):
  b=rd('low_hecke13_s2_completed.m'if q==13 and s==2 else f'low_hecke{q}_s{s}.m');assert not b[2];models.append(rd(f'low_hecke{q}_s{s}_model.m'))
  for r in b[0]:
   if r[0]in points:assert points[r[0]]==r
   else:points[r[0]]=r
  for r in b[1]:assert r[0]not in cases and r[4]==0;cases[r[0]]=r
 assert all(m==models[0]for m in models)
 wanted={r[0]:r for r in a[2]if(r[3]==0)==pure};assert set(wanted)==set(cases)
 for id,r in cases.items():assert r[:4]==wanted[id]
 targetrefs=sum(2+len(r[7])for r in points.values());data=[a[0],a[1],[points[i]for i in sorted(points)],[cases[i]for i in sorted(cases)]];(P/f'low_hecke{q}_merged.m').write_text(mag(data)+'\n');out.append({'hecke':q,'source_points':len(points),'source_cases':len(cases),'prime_references':targetrefs,'all_h0_zero':True})
(P/('low_hecke_merge_'+('_'.join(map(str,qs)))+'.json')).write_text(json.dumps({'status':'PASS','groups':out},indent=2)+'\n');print('D31B_C48_LOW_HECKE_MERGE_COMPLETE',out)
