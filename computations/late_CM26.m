SetColumns(0);SetQuitOnError(true);SetNthreads(1);
load "work/deg10/inputs/dcm_list_all_min_orders_X1_10k.m";
print "MINIMUM_CM_DEGREES",dcm_list_all_min_orders_X1[26],dcm_list_all_min_orders_X1[38],dcm_list_all_min_orders_X1[50];quit;
