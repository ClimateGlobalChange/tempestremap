#!/bin/sh

./GenerateTestData --mesh outCSne15.g --gll --np 2 --test 1 --out testdata_CSne15_np2_1.nc
./GenerateTestData --mesh outCSne15.g --gll --np 3 --test 1 --out testdata_CSne15_np3_1.nc
./GenerateTestData --mesh outCSne15.g --gll --np 4 --test 1 --out testdata_CSne15_np4_1.nc

./GenerateTestData --mesh outCSne15.g --gll --np 2 --test 2 --out testdata_CSne15_np2_2.nc
./GenerateTestData --mesh outCSne15.g --gll --np 3 --test 2 --out testdata_CSne15_np3_2.nc
./GenerateTestData --mesh outCSne15.g --gll --np 4 --test 2 --out testdata_CSne15_np4_2.nc

./GenerateTestData --mesh outCSne15.g --gll --np 2 --test 3 --out testdata_CSne15_np2_3.nc
./GenerateTestData --mesh outCSne15.g --gll --np 3 --test 3 --out testdata_CSne15_np3_3.nc
./GenerateTestData --mesh outCSne15.g --gll --np 4 --test 3 --out testdata_CSne15_np4_3.nc

