#!/usr/bin/env python3
"""Independent Sage/M4RI ranks, without native launch or output mutation."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path('work/deg10/continuation_20261030c/C65');t=time.monotonic();rr=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read());assert len(rr)==240
for c in rr:
 rows=[int(v,16)for v in c['rows']];n=c['columns'];assert len(rows)==121 and n>=239 and all(v.bit_length()<=n for v in rows);assert matrix(GF(2),[[v>>j&1 for j in range(n)]for v in rows]).rank()==121
print('DECIC65_R1_M4RI_RANK_PASS cases240 rank121 elapsed',time.monotonic()-t)

(H/'independent_m4ri.json').write_text(json.dumps({'status':'PASS','matrices':240,'all_ranks':121})+'\n')
