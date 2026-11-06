from math import gcd, isqrt
from pathlib import Path
import json
out = {}
for D in [-119, -319]:
    forms=[]
    for a in range(1,isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D)%(4*a): continue
            c=(b*b-D)//(4*a)
            if c<a or gcd(gcd(a,b),c)!=1: continue
            if (abs(b)==a or a==c) and b<0: continue
            forms.append([a,b,c])
    assert len(forms)==10
    out[str(D)]={"primitive_reduced_positive_forms":forms,"class_number":len(forms)}
assert all(pow(2,h//2,45) not in (1,44) for h in (2,10))
out["cycle_congruences"]={str(h):pow(2,h//2,45) for h in (2,10)}
Path("work/deg10/continuation_20261031a/C90/parent_CM_classnumbers.json").write_text(json.dumps(out,indent=2)+"\n")
print("C90_PARENT_CM_CLASSNUMBERS_COMPLETE",{D:out[str(D)]["class_number"] for D in [-119,-319]})
