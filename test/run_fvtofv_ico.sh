#!/bin/sh

time ../bin/GenerateTestData --mesh outCSne15.g --test 1 --out testdata_CSne15_1.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 2 --out testdata_CSne15_2.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 3 --out testdata_CSne15_3.nc

time ../bin/GenerateTestData --mesh outCSne30.g --test 1 --out testdata_CSne30_1.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 2 --out testdata_CSne30_2.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 3 --out testdata_CSne30_3.nc

time ../bin/GenerateTestData --mesh outCSne60.g --test 1 --out testdata_CSne60_1.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 2 --out testdata_CSne60_2.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 3 --out testdata_CSne60_3.nc

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_ICO72_np1_1.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_ICO72_np2_1.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_ICO72_np3_1.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_ICO72_np4_1.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_ICO72_np1_2.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_ICO72_np2_2.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_ICO72_np3_2.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_ICO72_np4_2.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_ICO72_np1_3.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_ICO72_np2_3.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_ICO72_np3_3.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outICO72.g --ov_mesh overlap_CSne15_ICO72.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_ICO72_np4_3.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_ICO72_np1_1.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_ICO72_np2_1.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_ICO72_np3_1.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_ICO72_np4_1.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_ICO72_np1_2.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_ICO72_np2_2.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_ICO72_np3_2.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_ICO72_np4_2.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_ICO72_np1_3.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_ICO72_np2_3.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_ICO72_np3_3.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO72.g --ov_mesh overlap_CSne30_ICO72.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_ICO72_np4_3.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_ICO72_np1_1.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_ICO72_np2_1.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_ICO72_np3_1.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_ICO72_np4_1.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_ICO72_np1_2.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_ICO72_np2_2.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_ICO72_np3_2.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_ICO72_np4_2.nc --in_np 4 --in_type fv --out_type fv

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_ICO72_np1_3.nc --in_np 1 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_ICO72_np2_3.nc --in_np 2 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_ICO72_np3_3.nc --in_np 3 --in_type fv --out_type fv
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO72.g --ov_mesh overlap_CSne60_ICO72.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_ICO72_np4_3.nc --in_np 4 --in_type fv --out_type fv

