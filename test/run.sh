#!/bin/sh
TEMPESTREMAP=.
$TEMPESTREMAP/GenerateCSMesh --res 64 --alt --file gravitySam.000000.3d.cubedSphere.g
$TEMPESTREMAP/GenerateRLLMesh --lon 256 --lat 128 --file gravitySam.000000.3d.latLon.g
$TEMPESTREMAP/GenerateOverlapMesh --a gravitySam.000000.3d.cubedSphere.g --b gravitySam.000000.3d.latLon.g --out gravitySam.000000.3d.overlap.g
$TEMPESTREMAP/GenerateOfflineMap --in_mesh gravitySam.000000.3d.cubedSphere.g --out_mesh gravitySam.000000.3d.latLon.g --ov_mesh gravitySam.000000.3d.overlap.g --in_np 2 --out_map gravitySam.000000.3d.mapping.nc 
