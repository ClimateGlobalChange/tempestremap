#!/bin/sh

time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --out testdata_CSne30t_np2_1.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --out testdata_CSne30t_np3_1.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --out testdata_CSne30t_np4_1.nc --gllint --np 4

time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --out testdata_CSne30t_np2_2.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --out testdata_CSne30t_np3_2.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --out testdata_CSne30t_np4_2.nc --gllint --np 4

time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --out testdata_CSne30t_np2_3.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --out testdata_CSne30t_np3_3.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --out testdata_CSne30t_np4_3.nc --gllint --np 4

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np2_1.nc --out_data testdata_CSne15_CSne30t_np2_np2_1.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np3_1.nc --out_data testdata_CSne15_CSne30t_np3_np3_1.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np4_1.nc --out_data testdata_CSne15_CSne30t_np4_np4_1.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np2_1.nc --out_data testdata_CSne30_CSne30t_np2_np2_1.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np3_1.nc --out_data testdata_CSne30_CSne30t_np3_np3_1.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np4_1.nc --out_data testdata_CSne30_CSne30t_np4_np4_1.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np2_1.nc --out_data testdata_CSne60_CSne30t_np2_np2_1.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np3_1.nc --out_data testdata_CSne60_CSne30t_np3_np3_1.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np4_1.nc --out_data testdata_CSne60_CSne30t_np4_np4_1.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble


time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np2_2.nc --out_data testdata_CSne15_CSne30t_np2_np2_2.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np3_2.nc --out_data testdata_CSne15_CSne30t_np3_np3_2.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np4_2.nc --out_data testdata_CSne15_CSne30t_np4_np4_2.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np2_2.nc --out_data testdata_CSne30_CSne30t_np2_np2_2.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np3_2.nc --out_data testdata_CSne30_CSne30t_np3_np3_2.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np4_2.nc --out_data testdata_CSne30_CSne30t_np4_np4_2.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np2_2.nc --out_data testdata_CSne60_CSne30t_np2_np2_2.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np3_2.nc --out_data testdata_CSne60_CSne30t_np3_np3_2.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np4_2.nc --out_data testdata_CSne60_CSne30t_np4_np4_2.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble


time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np2_3.nc --out_data testdata_CSne15_CSne30t_np2_np2_3.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np3_3.nc --out_data testdata_CSne15_CSne30t_np3_np3_3.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_np4_3.nc --out_data testdata_CSne15_CSne30t_np4_np4_3.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np2_3.nc --out_data testdata_CSne30_CSne30t_np2_np2_3.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np3_3.nc --out_data testdata_CSne30_CSne30t_np3_np3_3.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_np4_3.nc --out_data testdata_CSne30_CSne30t_np4_np4_3.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np2_3.nc --out_data testdata_CSne60_CSne30t_np2_np2_3.nc --in_np 2 --out_np 2 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np3_3.nc --out_data testdata_CSne60_CSne30t_np3_np3_3.nc --in_np 3 --out_np 3 --in_type cgll --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_np4_3.nc --out_data testdata_CSne60_CSne30t_np4_np4_3.nc --in_np 4 --out_np 4 --in_type cgll --out_type cgll --out_double --bubble

