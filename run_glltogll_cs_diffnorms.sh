#!/bin/sh

rm -rf testdata_glltogll_cs_diffnorms_1.txt
rm -rf testdata_glltogll_cs_diffnorms_2.txt
rm -rf testdata_glltogll_cs_diffnorms_3.txt

time ./GenerateTestData --mesh outCSne10t.g --test 1 --out testdata_CSne10t_np2_1.nc --gll --np 2
time ./GenerateTestData --mesh outCSne10t.g --test 1 --out testdata_CSne10t_np3_1.nc --gll --np 3
time ./GenerateTestData --mesh outCSne10t.g --test 1 --out testdata_CSne10t_np4_1.nc --gll --np 4

time ./GenerateTestData --mesh outCSne10t.g --test 2 --out testdata_CSne10t_np2_2.nc --gll --np 2
time ./GenerateTestData --mesh outCSne10t.g --test 2 --out testdata_CSne10t_np3_2.nc --gll --np 3
time ./GenerateTestData --mesh outCSne10t.g --test 2 --out testdata_CSne10t_np4_2.nc --gll --np 4

time ./GenerateTestData --mesh outCSne10t.g --test 3 --out testdata_CSne10t_np2_3.nc --gll --np 2
time ./GenerateTestData --mesh outCSne10t.g --test 3 --out testdata_CSne10t_np3_3.nc --gll --np 3
time ./GenerateTestData --mesh outCSne10t.g --test 3 --out testdata_CSne10t_np4_3.nc --gll --np 4

./CalculateDiffNorms --a testdata_CSne15_CSne10t_np2_np2_1.nc --b testdata_CSne10t_np2_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np3_np3_1.nc --b testdata_CSne10t_np3_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np4_np4_1.nc --b testdata_CSne10t_np4_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne30_CSne10t_np2_np2_1.nc --b testdata_CSne10t_np2_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np3_np3_1.nc --b testdata_CSne10t_np3_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np4_np4_1.nc --b testdata_CSne10t_np4_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt

./CalculateDiffNorms --a testdata_CSne60_CSne10t_np2_np2_1.nc --b testdata_CSne10t_np2_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np3_np3_1.nc --b testdata_CSne10t_np3_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np4_np4_1.nc --b testdata_CSne10t_np4_1.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_1.txt


./CalculateDiffNorms --a testdata_CSne15_CSne10t_np2_np2_2.nc --b testdata_CSne10t_np2_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np3_np3_2.nc --b testdata_CSne10t_np3_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np4_np4_2.nc --b testdata_CSne10t_np4_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne30_CSne10t_np2_np2_2.nc --b testdata_CSne10t_np2_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np3_np3_2.nc --b testdata_CSne10t_np3_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np4_np4_2.nc --b testdata_CSne10t_np4_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt

./CalculateDiffNorms --a testdata_CSne60_CSne10t_np2_np2_2.nc --b testdata_CSne10t_np2_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np3_np3_2.nc --b testdata_CSne10t_np3_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np4_np4_2.nc --b testdata_CSne10t_np4_2.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_2.txt


./CalculateDiffNorms --a testdata_CSne15_CSne10t_np2_np2_3.nc --b testdata_CSne10t_np2_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np3_np3_3.nc --b testdata_CSne10t_np3_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne15_CSne10t_np4_np4_3.nc --b testdata_CSne10t_np4_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne30_CSne10t_np2_np2_3.nc --b testdata_CSne10t_np2_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np3_np3_3.nc --b testdata_CSne10t_np3_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne30_CSne10t_np4_np4_3.nc --b testdata_CSne10t_np4_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt

./CalculateDiffNorms --a testdata_CSne60_CSne10t_np2_np2_3.nc --b testdata_CSne10t_np2_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np3_np3_3.nc --b testdata_CSne10t_np3_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt
./CalculateDiffNorms --a testdata_CSne60_CSne10t_np4_np4_3.nc --b testdata_CSne10t_np4_3.nc --mesh outCSne10t.g --outfile testdata_glltogll_cs_diffnorms_3.txt

