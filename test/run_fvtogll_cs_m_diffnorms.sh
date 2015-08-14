#!/bin/sh

rm -rf testdata_fvtogll_cs_mono_diffnorms_1.txt
rm -rf testdata_fvtogll_cs_mono_diffnorms_2.txt
rm -rf testdata_fvtogll_cs_mono_diffnorms_3.txt

time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --np 2 --out testdata_CSne30t_np2_1.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --np 3 --out testdata_CSne30t_np3_1.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 1 --np 4 --out testdata_CSne30t_np4_1.nc --gll

time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --np 2 --out testdata_CSne30t_np2_2.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --np 3 --out testdata_CSne30t_np3_2.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 2 --np 4 --out testdata_CSne30t_np4_2.nc --gll

time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --np 2 --out testdata_CSne30t_np2_3.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --np 3 --out testdata_CSne30t_np3_3.nc --gll
time ../bin/GenerateTestData --mesh outCSne30t.g --test 3 --np 4 --out testdata_CSne30t_np4_3.nc --gll

../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np2_1.nc --b testdata_CSne30t_np2_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np3_1.nc --b testdata_CSne30t_np3_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np4_1.nc --b testdata_CSne30t_np4_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np2_2.nc --b testdata_CSne30t_np2_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np3_2.nc --b testdata_CSne30t_np3_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np4_2.nc --b testdata_CSne30t_np4_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np2_3.nc --b testdata_CSne30t_np2_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np3_3.nc --b testdata_CSne30t_np3_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne15_CSne30t_mono_np4_3.nc --b testdata_CSne30t_np4_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np2_1.nc --b testdata_CSne30t_np2_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np3_1.nc --b testdata_CSne30t_np3_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np4_1.nc --b testdata_CSne30t_np4_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np2_2.nc --b testdata_CSne30t_np2_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np3_2.nc --b testdata_CSne30t_np3_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np4_2.nc --b testdata_CSne30t_np4_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np2_3.nc --b testdata_CSne30t_np2_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np3_3.nc --b testdata_CSne30t_np3_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne30_CSne30t_mono_np4_3.nc --b testdata_CSne30t_np4_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np2_1.nc --b testdata_CSne30t_np2_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np3_1.nc --b testdata_CSne30t_np3_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np4_1.nc --b testdata_CSne30t_np4_1.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_1.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np2_2.nc --b testdata_CSne30t_np2_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np3_2.nc --b testdata_CSne30t_np3_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np4_2.nc --b testdata_CSne30t_np4_2.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_2.txt --gll --np 4

../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np2_3.nc --b testdata_CSne30t_np2_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 2
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np3_3.nc --b testdata_CSne30t_np3_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 3
../bin/CalculateDiffNorms --a testdata_CSne60_CSne30t_mono_np4_3.nc --b testdata_CSne30t_np4_3.nc --mesh outCSne30t.g --outfile testdata_fvtogll_cs_mono_diffnorms_3.txt --gll --np 4

