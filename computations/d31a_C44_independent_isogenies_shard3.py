from pathlib import Path
import hashlib
SHARD=3
p=Path(__file__).with_name("d31a_C44_independent_isogenies_sharded_common.py")
s=p.read_text();assert hashlib.sha256(s.encode()).hexdigest()=="fa2bb7002c61088897460db0bb1b7be3d6c288b7856b11399f8986a96c1c8aca"
exec(compile(s,str(p),"exec"))
