from sage.all import GammaH,ModularSymbols
G=GammaH(68,[35,-1]);M=ModularSymbols(G,2,sign=0).cuspidal_subspace();W=M.atkin_lehner_operator(68).matrix()
for a in [1,3,5,7,9,11,13,15]:
 if a%17==0:continue
 D=M.diamond_bracket_matrix(a);A=W*D;print(a,A*A==1,int(A.trace()),'fixed',int(2-A.trace()),flush=True)
