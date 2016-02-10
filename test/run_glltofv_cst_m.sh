#!/bin/sh

time ../bin/GenerateOverlapMesh --a outCSne30t.g --b outRLL0.5deg.g --out overlap_CSne30t_RLL0.5deg.g

# CSne15 Source
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np2_1.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np2_2.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np2_3.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np3_1.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np3_2.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np3_3.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np4_1.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np4_2.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne15_CSne30t_mono_np4_3.nc --out_data testdata_CSne15_CSne30t_RLL0.5deg_mono_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double --mono

# CSne30 Source
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np2_1.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np2_2.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np2_3.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np3_1.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np3_2.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np3_3.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np4_1.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np4_2.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne30_CSne30t_mono_np4_3.nc --out_data testdata_CSne30_CSne30t_RLL0.5deg_mono_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double --mono

# CSne60 Source
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np2_1.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np2_2.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np2_3.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np3_1.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np3_2.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np3_3.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double --mono

time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np4_1.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np4_2.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double --mono
time ../bin/GenerateOfflineMap --in_mesh outCSne30t.g --out_mesh outRLL0.5deg.g --ov_mesh overlap_CSne30t_RLL0.5deg.g --var Psi --in_data testdata_CSne60_CSne30t_mono_np4_3.nc --out_data testdata_CSne60_CSne30t_RLL0.5deg_mono_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double --mono





