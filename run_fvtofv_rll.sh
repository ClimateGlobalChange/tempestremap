#!/bin/sh

time ./GenerateTestData --mesh outCSne15.g --test 1 --out testdata_CSne15_1.nc
time ./GenerateTestData --mesh outCSne15.g --test 2 --out testdata_CSne15_2.nc
time ./GenerateTestData --mesh outCSne15.g --test 3 --out testdata_CSne15_3.nc

time ./GenerateTestData --mesh outCSne30.g --test 1 --out testdata_CSne30_1.nc
time ./GenerateTestData --mesh outCSne30.g --test 2 --out testdata_CSne30_2.nc
time ./GenerateTestData --mesh outCSne30.g --test 3 --out testdata_CSne30_3.nc

time ./GenerateTestData --mesh outCSne60.g --test 1 --out testdata_CSne60_1.nc
time ./GenerateTestData --mesh outCSne60.g --test 2 --out testdata_CSne60_2.nc
time ./GenerateTestData --mesh outCSne60.g --test 3 --out testdata_CSne60_3.nc

time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_RLL2deg_np1_1.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_RLL2deg_np2_1.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_RLL2deg_np3_1.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_RLL2deg_np4_1.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_RLL2deg_np1_2.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_RLL2deg_np2_2.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_RLL2deg_np3_2.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_RLL2deg_np4_2.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_RLL2deg_np1_3.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_RLL2deg_np2_3.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_RLL2deg_np3_3.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne15_RLL2deg.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_RLL2deg_np4_3.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_RLL2deg_np1_1.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_RLL2deg_np2_1.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_RLL2deg_np3_1.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_RLL2deg_np4_1.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_RLL2deg_np1_2.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_RLL2deg_np2_2.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_RLL2deg_np3_2.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_RLL2deg_np4_2.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_RLL2deg_np1_3.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_RLL2deg_np2_3.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_RLL2deg_np3_3.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne30_RLL2deg.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_RLL2deg_np4_3.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_RLL2deg_np1_1.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_RLL2deg_np2_1.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_RLL2deg_np3_1.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_RLL2deg_np4_1.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_RLL2deg_np1_2.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_RLL2deg_np2_2.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_RLL2deg_np3_2.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_RLL2deg_np4_2.nc --np 4

time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_RLL2deg_np1_3.nc --np 1
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_RLL2deg_np2_3.nc --np 2
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_RLL2deg_np3_3.nc --np 3
time ./GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outRLL2deg.g --ov_mesh overlap_CSne60_RLL2deg.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_RLL2deg_np4_3.nc --np 4

