from pathlib import Path
p=Path("work/deg10/continuation_20261031a/C46/independent_quotient_arithmetic.py")
exec(compile(p.read_text(),str(p),"exec"),{"__file__":str(p.resolve()),"__name__":"__main__"})
