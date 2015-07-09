#!/bin/sh

rm -rf testdata_fvtofv_rll_diffnorms_1.txt
rm -rf testdata_fvtofv_rll_diffnorms_2.txt
rm -rf testdata_fvtofv_rll_diffnorms_3.txt

time ../bin/GenerateTestData --mesh outRLL2deg.g --test 1 --out testdata_RLL2deg_1.nc
time ../bin/GenerateTestData --mesh outRLL2deg.g --test 2 --out testdata_RLL2deg_2.nc
time ../bin/GenerateTestData --mesh outRLL2deg.g --test 3 --out testdata_RLL2deg_3.nc

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np1_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np1_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np1_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np1_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np1_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np1_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np1_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np1_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np1_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_fvtofv_rll_diffnorms_3.txt

