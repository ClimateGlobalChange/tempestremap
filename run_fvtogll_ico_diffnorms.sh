#!/bin/sh

rm -rf testdata_fvtogll_ico_diffnorms_1.txt
rm -rf testdata_fvtogll_ico_diffnorms_2.txt
rm -rf testdata_fvtogll_ico_diffnorms_3.txt

time ./GenerateTestData --mesh outCSne30.g --test 1 --out testdata_CSne30_np4_1.nc --gll
time ./GenerateTestData --mesh outCSne30.g --test 2 --out testdata_CSne30_np4_2.nc --gll
time ./GenerateTestData --mesh outCSne30.g --test 3 --out testdata_CSne30_np4_3.nc --gll

./CalculateDiffNorms --a testdata_ICO18_CSne30_np1_1.nc --b testdata_CSne30_np4_1.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_1.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np2_1.nc --b testdata_CSne30_np4_1.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_1.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np3_1.nc --b testdata_CSne30_np4_1.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_1.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np4_1.nc --b testdata_CSne30_np4_1.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_1.txt --gll

./CalculateDiffNorms --a testdata_ICO18_CSne30_np1_2.nc --b testdata_CSne30_np4_2.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_2.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np2_2.nc --b testdata_CSne30_np4_2.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_2.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np3_2.nc --b testdata_CSne30_np4_2.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_2.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np4_2.nc --b testdata_CSne30_np4_2.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_2.txt --gll

./CalculateDiffNorms --a testdata_ICO18_CSne30_np1_3.nc --b testdata_CSne30_np4_3.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_3.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np2_3.nc --b testdata_CSne30_np4_3.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_3.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np3_3.nc --b testdata_CSne30_np4_3.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_3.txt --gll
./CalculateDiffNorms --a testdata_ICO18_CSne30_np4_3.nc --b testdata_CSne30_np4_3.nc --mesh outCSne30.g --outfile testdata_fvtogll_ico_diffnorms_3.txt --gll

