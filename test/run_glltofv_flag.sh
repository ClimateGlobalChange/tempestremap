#!/bin/sh

time ../bin/GenerateOverlapMesh --a outFLAG16.g --b outRLL1deg.g --out overlap_FLAG16_RLL1deg.g

# CSne15 Source
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np2_np2_1.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np2_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np2_np2_2.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np2_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np2_np2_3.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np2_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np3_np3_1.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np3_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np3_np3_2.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np3_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np3_np3_3.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np3_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np4_np4_1.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np4_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np4_np4_2.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np4_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne15_FLAG16_np4_np4_3.nc --out_data testdata_CSne15_FLAG16_RLL1deg_np4_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double

# CSne30 Source
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np2_np2_1.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np2_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np2_np2_2.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np2_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np2_np2_3.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np2_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np3_np3_1.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np3_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np3_np3_2.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np3_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np3_np3_3.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np3_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np4_np4_1.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np4_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np4_np4_2.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np4_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne30_FLAG16_np4_np4_3.nc --out_data testdata_CSne30_FLAG16_RLL1deg_np4_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double

# CSne60 Source
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np2_np2_1.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np2_np2_1.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np2_np2_2.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np2_np2_2.nc --in_type cgll --in_np 2 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np2_np2_3.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np2_np2_3.nc --in_type cgll --in_np 2 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np3_np3_1.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np3_np3_1.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np3_np3_2.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np3_np3_2.nc --in_type cgll --in_np 3 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np3_np3_3.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np3_np3_3.nc --in_type cgll --in_np 3 --bubble --out_double

time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np4_np4_1.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np4_np4_1.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np4_np4_2.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np4_np4_2.nc --in_type cgll --in_np 4 --bubble --out_double
time ../bin/GenerateOfflineMap --in_mesh outFLAG16.g --out_mesh outRLL1deg.g --ov_mesh overlap_FLAG16_RLL1deg.g --var Psi --in_data testdata_CSne60_FLAG16_np4_np4_3.nc --out_data testdata_CSne60_FLAG16_RLL1deg_np4_np4_3.nc --in_type cgll --in_np 4 --bubble --out_double







