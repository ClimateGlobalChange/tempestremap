#!/bin/sh

rm -rf testdata_glltofv_rll_mono_diffnorms_1.txt
rm -rf testdata_glltofv_rll_mono_diffnorms_2.txt
rm -rf testdata_glltofv_rll_mono_diffnorms_3.txt

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np2_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np3_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np4_1.nc --b testdata_RLL2deg_1.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_1.txt

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np2_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np3_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np4_2.nc --b testdata_RLL2deg_2.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_2.txt

../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne15_RLL2deg_mono_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt

../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne30_RLL2deg_mono_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt

../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np2_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np3_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt
../bin/CalculateDiffNorms --a testdata_CSne60_RLL2deg_mono_np4_3.nc --b testdata_RLL2deg_3.nc --mesh outRLL2deg.g --outfile testdata_glltofv_rll_mono_diffnorms_3.txt

