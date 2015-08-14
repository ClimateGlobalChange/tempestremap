#!/bin/sh

time ../bin/GenerateOverlapMesh --a outCSne15.g --b outCSne30t.g --out overlap_CSne15_CSne30t.g
time ../bin/GenerateOverlapMesh --a outCSne30.g --b outCSne30t.g --out overlap_CSne30_CSne30t.g
time ../bin/GenerateOverlapMesh --a outCSne60.g --b outCSne30t.g --out overlap_CSne60_CSne30t.g

time ../bin/GenerateTestData --mesh outCSne15.g --test 1 --out testdata_CSne15_1.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 2 --out testdata_CSne15_2.nc
time ../bin/GenerateTestData --mesh outCSne15.g --test 3 --out testdata_CSne15_3.nc

time ../bin/GenerateTestData --mesh outCSne30.g --test 1 --out testdata_CSne30_1.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 2 --out testdata_CSne30_2.nc
time ../bin/GenerateTestData --mesh outCSne30.g --test 3 --out testdata_CSne30_3.nc

time ../bin/GenerateTestData --mesh outCSne60.g --test 1 --out testdata_CSne60_1.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 2 --out testdata_CSne60_2.nc
time ../bin/GenerateTestData --mesh outCSne60.g --test 3 --out testdata_CSne60_3.nc

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_CSne30t_np1_1.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_CSne30t_np2_1.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_CSne30t_np3_1.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_1.nc --out_data testdata_CSne15_CSne30t_np4_1.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_CSne30t_np1_2.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_CSne30t_np2_2.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_CSne30t_np3_2.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_2.nc --out_data testdata_CSne15_CSne30t_np4_2.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_CSne30t_np1_3.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_CSne30t_np2_3.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_CSne30t_np3_3.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne15.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne15_CSne30t.g --var Psi --in_data testdata_CSne15_3.nc --out_data testdata_CSne15_CSne30t_np4_3.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_CSne30t_np1_1.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_CSne30t_np2_1.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_CSne30t_np3_1.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_1.nc --out_data testdata_CSne30_CSne30t_np4_1.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_CSne30t_np1_2.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_CSne30t_np2_2.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_CSne30t_np3_2.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_2.nc --out_data testdata_CSne30_CSne30t_np4_2.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_CSne30t_np1_3.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_CSne30t_np2_3.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_CSne30t_np3_3.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne30_CSne30t.g --var Psi --in_data testdata_CSne30_3.nc --out_data testdata_CSne30_CSne30t_np4_3.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_CSne30t_np1_1.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_CSne30t_np2_1.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_CSne30t_np3_1.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_1.nc --out_data testdata_CSne60_CSne30t_np4_1.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_CSne30t_np1_2.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_CSne30t_np2_2.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_CSne30t_np3_2.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_2.nc --out_data testdata_CSne60_CSne30t_np4_2.nc --in_np 4 --out_type cgll --out_double --bubble

time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_CSne30t_np1_3.nc --in_np 1 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_CSne30t_np2_3.nc --in_np 2 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_CSne30t_np3_3.nc --in_np 3 --out_type cgll --out_double --bubble
time ../bin/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outCSne30t.g --ov_mesh overlap_CSne60_CSne30t.g --var Psi --in_data testdata_CSne60_3.nc --out_data testdata_CSne60_CSne30t_np4_3.nc --in_np 4 --out_type cgll --out_double --bubble

