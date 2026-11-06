// Recompute the degree9 certificate without overwriting its frozen evidence.
JetFile:="work/degree9/realizations/continuation/C63/reproduced_cartier_data.m";
fresh:=Open(JetFile,"w");delete fresh;
load "work/degree9/realizations/continuation/C63/cartier_rr_D9.m";
