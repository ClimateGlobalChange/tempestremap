#!/bin/sh

time ../bin/GenerateTestData --mesh outCSne15.g --test 1 --gll --np 2 --out testdata_CSne15_np2_1.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 1 --gll --np 3 --out testdata_CSne15_np3_1.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 1 --gll --np 4 --out testdata_CSne15_np4_1.nc

time ../bin/GenerateTestData --mesh outCSne15.g --test 2 --gll --np 2 --out testdata_CSne15_np2_2.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 2 --gll --np 3 --out testdata_CSne15_np3_2.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 2 --gll --np 4 --out testdata_CSne15_np4_2.nc

time ../bin/GenerateTestData --mesh outCSne15.g --test 3 --gll --np 2 --out testdata_CSne15_np2_3.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 3 --gll --np 3 --out testdata_CSne15_np3_3.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 3 --gll --np 4 --out testdata_CSne15_np4_3.nc

time ../bin/GenerateTestData --mesh outCSne30.g --test 1 --gll --np 2 --out testdata_CSne30_np2_1.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 1 --gll --np 3 --out testdata_CSne30_np3_1.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 1 --gll --np 4 --out testdata_CSne30_np4_1.nc

time ../bin/GenerateTestData --mesh outCSne30.g --test 2 --gll --np 2 --out testdata_CSne30_np2_2.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 2 --gll --np 3 --out testdata_CSne30_np3_2.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 2 --gll --np 4 --out testdata_CSne30_np4_2.nc

time ../bin/GenerateTestData --mesh outCSne30.g --test 3 --gll --np 2 --out testdata_CSne30_np2_3.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 3 --gll --np 3 --out testdata_CSne30_np3_3.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 3 --gll --np 4 --out testdata_CSne30_np4_3.nc

time ../bin/GenerateTestData --mesh outCSne60.g --test 1 --gll --np 2 --out testdata_CSne60_np2_1.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 1 --gll --np 3 --out testdata_CSne60_np3_1.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 1 --gll --np 4 --out testdata_CSne60_np4_1.nc

time ../bin/GenerateTestData --mesh outCSne60.g --test 2 --gll --np 2 --out testdata_CSne60_np2_2.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 2 --gll --np 3 --out testdata_CSne60_np3_2.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 2 --gll --np 4 --out testdata_CSne60_np4_2.nc

time ../bin/GenerateTestData --mesh outCSne60.g --test 3 --gll --np 2 --out testdata_CSne60_np2_3.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 3 --gll --np 3 --out testdata_CSne60_np3_3.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 3 --gll --np 4 --out testdata_CSne60_np4_3.nc

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np2_1.nc --out_data testdata_CSne15_RLL1deg_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np3_1.nc --out_data testdata_CSne15_RLL1deg_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np4_1.nc --out_data testdata_CSne15_RLL1deg_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np2_2.nc --out_data testdata_CSne15_RLL1deg_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np3_2.nc --out_data testdata_CSne15_RLL1deg_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np4_2.nc --out_data testdata_CSne15_RLL1deg_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np2_3.nc --out_data testdata_CSne15_RLL1deg_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np3_3.nc --out_data testdata_CSne15_RLL1deg_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne15_RLL1deg.g --var Psi --in_data testdata_CSne15_np4_3.nc --out_data testdata_CSne15_RLL1deg_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np2_1.nc --out_data testdata_CSne30_RLL1deg_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np3_1.nc --out_data testdata_CSne30_RLL1deg_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np4_1.nc --out_data testdata_CSne30_RLL1deg_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np2_2.nc --out_data testdata_CSne30_RLL1deg_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np3_2.nc --out_data testdata_CSne30_RLL1deg_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np4_2.nc --out_data testdata_CSne30_RLL1deg_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np2_3.nc --out_data testdata_CSne30_RLL1deg_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np3_3.nc --out_data testdata_CSne30_RLL1deg_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne30_RLL1deg.g --var Psi --in_data testdata_CSne30_np4_3.nc --out_data testdata_CSne30_RLL1deg_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np2_1.nc --out_data testdata_CSne60_RLL1deg_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np3_1.nc --out_data testdata_CSne60_RLL1deg_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np4_1.nc --out_data testdata_CSne60_RLL1deg_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np2_2.nc --out_data testdata_CSne60_RLL1deg_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np3_2.nc --out_data testdata_CSne60_RLL1deg_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np4_2.nc --out_data testdata_CSne60_RLL1deg_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np2_3.nc --out_data testdata_CSne60_RLL1deg_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np3_3.nc --out_data testdata_CSne60_RLL1deg_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL1deg.g --ov_mesh overlap_CSne60_RLL1deg.g --var Psi --in_data testdata_CSne60_np4_3.nc --out_data testdata_CSne60_RLL1deg_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double

