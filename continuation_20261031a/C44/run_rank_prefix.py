from pathlib import Path
import re,subprocess,json,hashlib,time
H=Path(__file__).resolve().parent;R=H.parents[3]
log=R/'work/deg10/jobs/d31a_C44_cusp_jets.log'
assert 'D31A_C44_CUSP_JETS_COMPLETE places=24' in log.read_text()
s=(H/'cusp_jets.m').read_text();ds=list(map(int,re.search(r'Ds:=\[(.*?)\];',s,re.S).group(1).replace(',',' ').split()))
bs=[]
for x in re.finditer(r'Append\(~BasisJets,Matrix\(GF\(3\),36,(\d+),\[(.*?)\]\)\);',s,re.S):
 d=int(x.group(1));v=list(map(int,re.findall(r'\b[012]\b',x.group(2))));assert len(v)==36*d
 bs.append([v[i*d:(i+1)*d]for i in range(36)])
assert len(ds)==len(bs)==24
q=H/'cusp_columns.txt'
with q.open('w')as f:
 print(24,36,10,file=f)
 for d,b in zip(ds,bs):
  nc=len(b[0]);assert nc==d*(10//d);print(d,nc,file=f)
  for j in range(nc):print(sum((b[i][j]==1)<<i for i in range(36)),sum((b[i][j]==2)<<i for i in range(36)),file=f)
subprocess.run(['g++','-O3','-std=c++17',str(H/'rank_prefix.cpp'),'-o',str(H/'rank_prefix')],check=True)
t=time.time();p=subprocess.run([str(H/'rank_prefix'),str(q),str(H/'cusp_rank_exceptions.jsonl')],check=True,capture_output=True,text=True,timeout=300);print(p.stdout.strip())
(H/'cusp_prefix_replay.json').write_text(json.dumps(dict(status='PASS',summary=p.stdout.strip(),seconds=time.time()-t,
 source_sha256=hashlib.sha256(s.encode()).hexdigest(),cpp_sha256=hashlib.sha256((H/'rank_prefix.cpp').read_bytes()).hexdigest()),indent=2)+chr(10))
