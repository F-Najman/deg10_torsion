from sage.all import *
from pathlib import Path
import ast,re,json,hashlib
ROOT=Path('work/deg10/continuation_20261030c/C91_CM_control')
text=(ROOT/'inputs/records91_8.m').read_text()
def lit(name):return ast.literal_eval(re.search(name+r':=(.*?);',text,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
R=PolynomialRing(GF(2),'T');pol=lit('FieldPolynomial');K=GF(256,'a',modulus=R(pol));a=K.gen();P=PolynomialRing(K,'x');x=P.gen();out=[]
for i,rec in enumerate(lit('Records'),1):
 E=EllipticCurve(K,[R(v)(a)for v in rec[1]]);Q=E(R(rec[2][0])(a),R(rec[2][1])(a));assert E.cardinality()==273 and Q.order()==91
 roots=(x*x+x+1).roots(multiplicities=False);zz=[u for u in roots if E(u*Q[0],Q[1])==16*Q];assert len(zz)==1;u=zz[0]
 Q7=13*Q;S13=7*Q;assert Q7.order()==7 and S13.order()==13 and E(u*Q7[0],Q7[1])==2*Q7 and E(u*S13[0],S13[1])==3*S13
 ker=prod(x-(j*S13)[0] for j in range(1,7));assert prod(x-u*(j*S13)[0] for j in range(1,7))==ker
 out.append({'source_twist':i,'omega_scalar_on_P91':16,'omega_scalar_on_P7':2,'omega_scalar_on_S13':3,'cyclic13_kernel_invariant':True})
(ROOT/'independent_CM_level7_binding.json').write_text(json.dumps({'status':'PASS','twists':2,'geometric_tau2_fix_of_all_point7_cyclic13_images':True,'records':out,'input_sha256':hashlib.sha256((ROOT/'inputs/records91_8.m').read_bytes()).hexdigest()},indent=2)+'\n')
print('C91_INDEPENDENT_SAGE_CM_POINT7_KERNEL13_BINDING_PASS',2)
