load "work/deg10/continuation_20261031a/C62/class_common_v2.m";
qq:=Places(F,10);load "work/deg10/continuation_20261031a/C62/degree10_classes.m";
a:=ClassTests[1];cl,wit:=SmoothClass(qq[a[1]]);assert [Eltseq(v):v in cl] eq a[2];
print "D31A_C62_OPTIMIZED_SMOOTHING_VERIFIED",a[1],Cputime(tt);quit;
