from pathlib import Path
p=Path("work/deg10/continuation_20261031a/C44/run_rank_prefix.py")
exec(compile(p.read_text(),str(p),"exec"),{"__file__":str(p.resolve()),"__name__":"__main__"})
