from pathlib import Path
import ast,json,math
from sage.all import ZZ,QQ
P=Path('work/deg10/C226_resolution_20261030')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
primes=[1000000007,1000000009,1000000021,1000000033];data=[read(P/f'finite_units_{p}.m') for p in primes];out=[]
def lift(res):
 a=0;m=1
 for b,p in zip(res,primes):a+=((b-a)*pow(m,-1,p)%p)*m;m*=p
 q=ZZ(a).rational_reconstruction(ZZ(m));assert all((int(q.numerator())-int(q.denominator())*b)%p==0 for b,p in zip(res,primes));return [int(q.numerator()),int(q.denominator())]
for j in range(2):
 rows=[d[1][j] for d in data];idx=rows[0][0];assert all(r[0]==idx and r[3:]==rows[0][3:] and len(r[1])==len(rows[0][1]) for r in rows)
 den=[lift([r[1][i] for r in rows]) for i in range(len(rows[0][1]))];nums=[]
 for k in range(12):
  assert all(len(r[2][k])==len(rows[0][2][k]) for r in rows);nums.append([lift([r[2][k][i] for r in rows]) for i in range(len(rows[0][2][k]))])
 out.append([idx,den,nums]);co=[v for pol in [den]+nums for v in pol];print('Q226_RECONSTRUCTED_UNIT',idx,'max_num_bits',max(abs(v[0]).bit_length() for v in co),'max_den_bits',max(v[1].bit_length() for v in co),flush=True)
def enc(x):
 if isinstance(x,list):
  a,b=('[',']') if all(isinstance(v,int) for v in x) else ('[*','*]');return a+','.join(enc(v) for v in x)+b
 return str(x)
(P/'reconstructed_units.m').write_text(enc(out)+'\n');(P/'unit_reconstruction.json').write_text(json.dumps({'status':'CANDIDATE; EXACT Q VALIDATION REQUIRED','primes':primes,'classes':[r[0] for r in out]},indent=2)+'\n')
print('Q226_FUNCTION_RECONSTRUCTION_COMPLETE',flush=True)
