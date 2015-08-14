#!/bin/sh

time ./GenerateTestData --mesh outRLL1deg.g --test 1 --out testdata_RLL1deg_1.nc
time ./GenerateTestData --mesh outRLL1deg.g --test 2 --out testdata_RLL1deg_2.nc
time ./GenerateTestData --mesh outRLL1deg.g --test 3 --out testdata_RLL1deg_3.nc

time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_1.nc --out_data testdata_RLL1deg_CSne30_np1_1.nc --in_np 1 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_1.nc --out_data testdata_RLL1deg_CSne30_np2_1.nc --in_np 2 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_1.nc --out_data testdata_RLL1deg_CSne30_np3_1.nc --in_np 3 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_1.nc --out_data testdata_RLL1deg_CSne30_np4_1.nc --in_np 4 --out_type cgll --out_double --bubble

time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_2.nc --out_data testdata_RLL1deg_CSne30_np1_2.nc --in_np 1 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_2.nc --out_data testdata_RLL1deg_CSne30_np2_2.nc --in_np 2 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_2.nc --out_data testdata_RLL1deg_CSne30_np3_2.nc --in_np 3 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_2.nc --out_data testdata_RLL1deg_CSne30_np4_2.nc --in_np 4 --out_type cgll --out_double --bubble

time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_3.nc --out_data testdata_RLL1deg_CSne30_np1_3.nc --in_np 1 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_3.nc --out_data testdata_RLL1deg_CSne30_np2_3.nc --in_np 2 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_3.nc --out_data testdata_RLL1deg_CSne30_np3_3.nc --in_np 3 --out_type cgll --out_double --bubble
time ./GenerateOfflineMap --in_mesh outRLL1deg.g --out_mesh outCSne30.g --ov_mesh overlap_RLL1deg_CSne30.g --var Psi --in_data testdata_RLL1deg_3.nc --out_data testdata_RLL1deg_CSne30_np4_3.nc --in_np 4 --out_type cgll --out_double --bubble

