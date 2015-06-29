#!/bin/sh

rm -rf testdata_glltofv_rll_diffnorms_1.txt
rm -rf testdata_glltofv_rll_diffnorms_2.txt
rm -rf testdata_glltofv_rll_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne15_RLL55km_np2_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np3_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np4_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne30_RLL55km_np2_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np3_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np4_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne60_RLL55km_np2_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np3_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np4_1.nc --b testdata_RLL55km_1.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne15_RLL55km_np2_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np3_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np4_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne30_RLL55km_np2_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np3_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np4_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne60_RLL55km_np2_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np3_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np4_2.nc --b testdata_RLL55km_2.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne15_RLL55km_np2_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np3_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL55km_np4_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne30_RLL55km_np2_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np3_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL55km_np4_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne60_RLL55km_np2_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np3_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL55km_np4_3.nc --b testdata_RLL55km_3.nc --mesh outRLL55km.g --outfile testdata_glltofv_rll_diffnorms_3.txt

