#!/bin/sh

rm -rf testdata_glltogll_cs_diffnorms_1.txt
rm -rf testdata_glltogll_cs_diffnorms_2.txt
rm -rf testdata_glltogll_cs_diffnorms_3.txt

time ../bin/GenerateTestData --mesh outFLAG16.g --test 1 --out testdata_FLAG16_np2_1.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outFLAG16.g --test 1 --out testdata_FLAG16_np3_1.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outFLAG16.g --test 1 --out testdata_FLAG16_np4_1.nc --gllint --np 4

time ../bin/GenerateTestData --mesh outFLAG16.g --test 2 --out testdata_FLAG16_np2_2.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outFLAG16.g --test 2 --out testdata_FLAG16_np3_2.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outFLAG16.g --test 2 --out testdata_FLAG16_np4_2.nc --gllint --np 4

time ../bin/GenerateTestData --mesh outFLAG16.g --test 3 --out testdata_FLAG16_np2_3.nc --gllint --np 2
time ../bin/GenerateTestData --mesh outFLAG16.g --test 3 --out testdata_FLAG16_np3_3.nc --gllint --np 3
time ../bin/GenerateTestData --mesh outFLAG16.g --test 3 --out testdata_FLAG16_np4_3.nc --gllint --np 4

../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np2_np2_1.nc --b testdata_FLAG16_np2_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np3_np3_1.nc --b testdata_FLAG16_np3_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np4_np4_1.nc --b testdata_FLAG16_np4_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np2_np2_1.nc --b testdata_FLAG16_np2_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np3_np3_1.nc --b testdata_FLAG16_np3_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np4_np4_1.nc --b testdata_FLAG16_np4_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np2_np2_1.nc --b testdata_FLAG16_np2_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np3_np3_1.nc --b testdata_FLAG16_np3_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np4_np4_1.nc --b testdata_FLAG16_np4_1.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_1.txt --np 4 --gll


../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np2_np2_2.nc --b testdata_FLAG16_np2_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np3_np3_2.nc --b testdata_FLAG16_np3_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np4_np4_2.nc --b testdata_FLAG16_np4_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np2_np2_2.nc --b testdata_FLAG16_np2_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np3_np3_2.nc --b testdata_FLAG16_np3_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np4_np4_2.nc --b testdata_FLAG16_np4_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np2_np2_2.nc --b testdata_FLAG16_np2_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np3_np3_2.nc --b testdata_FLAG16_np3_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np4_np4_2.nc --b testdata_FLAG16_np4_2.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_2.txt --np 4 --gll


../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np2_np2_3.nc --b testdata_FLAG16_np2_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np3_np3_3.nc --b testdata_FLAG16_np3_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne15_FLAG16_np4_np4_3.nc --b testdata_FLAG16_np4_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np2_np2_3.nc --b testdata_FLAG16_np2_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np3_np3_3.nc --b testdata_FLAG16_np3_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne30_FLAG16_np4_np4_3.nc --b testdata_FLAG16_np4_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 4 --gll

../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np2_np2_3.nc --b testdata_FLAG16_np2_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 2 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np3_np3_3.nc --b testdata_FLAG16_np3_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 3 --gll
../bin/CalculateDiffNorms --a testdata_CSne60_FLAG16_np4_np4_3.nc --b testdata_FLAG16_np4_3.nc --mesh outFLAG16.g --outfile testdata_glltogll_cs_diffnorms_3.txt --np 4 --gll

