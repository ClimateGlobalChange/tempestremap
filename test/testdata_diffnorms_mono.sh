#!/bin/sh

rm -rf testdata_diffnorms_mono_1.txt
rm -rf testdata_diffnorms_mono_2.txt
rm -rf testdata_diffnorms_mono_3.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_mono_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_1.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_mono_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_2.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_mono_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_mono_3.txt

