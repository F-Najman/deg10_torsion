from pathlib import Path
import hashlib
SHARD=3
p=Path(__file__).with_name("d31a_C44_independent_isogenies_sharded_common_v2.py")
s=p.read_text();assert hashlib.sha256(s.encode()).hexdigest()=="cc840d4b273d75be5b4aa9df6525c76a3159780251db27912330d987b6fc7318"
exec(compile(s,str(p),"exec"))
