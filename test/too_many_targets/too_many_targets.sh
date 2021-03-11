#!/bin/bash

CSRES=12

rm *.g
rm map*.nc

TEMPESTREMAPDIR=../../bin

$TEMPESTREMAPDIR/GenerateRLLMesh --lon 360 --lat 180 --file outRLL1deg.g

$TEMPESTREMAPDIR/GenerateCSMesh --res ${CSRES} --file outCSne${CSRES}.g

$TEMPESTREMAPDIR/GenerateOverlapMesh --a outRLL1deg.g --b outCSne${CSRES}.g --out ov_RLL1deg_CSne${CSRES}.g

$TEMPESTREMAPDIR/GenerateOfflineMap --in_mesh outCSne${CSRES}.g --out_mesh outRLL1deg.g --ov_mesh ov_RLL1deg_CSne${CSRES}.g --in_type cgll --in_np 4 --out_type fv --out_map map_CSne${CSRES}np4_RLL1deg.nc
