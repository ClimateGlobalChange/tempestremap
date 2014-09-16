#!/bin/sh

rm -rf testdata_diffnorms_1.txt
rm -rf testdata_diffnorms_2.txt
rm -rf testdata_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_1.nc --b testdata_RLL1deg_1.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_2.nc --b testdata_RLL1deg_2.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np2_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np3_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_RLL1deg_np4_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np2_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np3_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_RLL1deg_np4_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np2_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np3_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_RLL1deg_np4_3.nc --b testdata_RLL1deg_3.nc --mesh outRLL1deg.g --outfile testdata_diffnorms_3.txt

