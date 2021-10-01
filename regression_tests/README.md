***
# Regression Testing for TempestRemap Library
## Overview:
- Paths used for regression are relative to the directory <yourpath>/tempestremap/regression_tests/
    - Python calls must be made from this directory
    - On a newly cloned repository:
        - Generate baseline using ```python regression_tests.py -g```  generates ```baseline_data.pkl``` and ```timing_data.pkl```
        - Now run the regression as required: ```python regression_tests.py ``` generates ```baseline_data_repeat.pkl``` and ```regtime.pkl```
    - There are 5 major commands in the regression workflow
      - Generate mesh (called twice)
        - GenerateICOMesh
        - GenerateCSMesh
        - GenerateRLLMesh
      - GenerateTestData
      - GenerateOfflineMap
      - GenerateOverlapMesh
      - ApplyOfflineMap
    - % difference b/w source and target is calculated and reported
    - The test_matrix.ini file  specifies the command line options for each of the 6 commands
      - test_matrix.ini
```
            id   srcmesh   tgtmesh   order correctareas monotone test
            0    cs-30     icod-60    1         0          0       1
            1    cs-11     icod-30    1         1          0       2
            2    cs-30     icod-60    2         0          0       1
            3    cs-30     icod-60    3         0          0       1
            4    cs-30     icod-60    4         0          0       1
            5    cs-30     icod-60    1         0          0       2
            6    cs-11     rll-30-60  1         1          1       1
``` 

*** 
## Directory Structure
- regression_tests
    - baseline: pickle files for baseline runs, baseline timing and regression time.
    - data: test data files
    - maps: map files
    - meshes: mesh files

### Pickle File Structures:
#### baseline/baseline_data.pkl
```
id     source     target         error
0   0  25.132741  25.132741  1.130864e-13
1   1  25.132733  25.132733  1.979012e-13
2   2  25.132741  25.132741  1.130864e-13
3   3  25.132741  25.132741  1.074321e-12
4   4  25.132741  25.132741  3.251234e-13
5   5  25.132741  25.132741  3.533950e-13
6   6  25.132741  25.132741  4.240740e-14
```

#### baseline/timing_data.pkl
```
                                              command  average_time  num_runs
0   GenerateICOMesh --dual --res 60 --file meshes/...      0.658594         6
1   GenerateCSMesh --file meshes/outCSMesh-11.g --...      0.128430         6
2   GenerateCSMesh --file meshes/outCSMesh-30.g --...      0.083597         6
3   GenerateICOMesh --dual --res 30 --file meshes/...      0.234755         6
4   GenerateRLLMesh --lon 30 --lat 60 --file meshe...      0.099860         6
5   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.399087         6
6   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.415866         6
7   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.411952         6
8   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.445240         6
9   GenerateTestData --mesh meshes/outCSMesh-11.g ...      0.147044         6
10  GenerateTestData --mesh meshes/outCSMesh-11.g ...      0.118803         6
11  GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.395406         6
12  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...      0.493001         6
13  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...      1.508532         6
14  GenerateOverlapMesh --a meshes/outCSMesh-30.g ...      5.214701         6
15  GenerateOfflineMap --in_mesh meshes/outCSMesh-...     10.623075         6
16  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      7.118552         6
17  GenerateOfflineMap --in_mesh meshes/outCSMesh-...     16.365483         6
18  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      0.392716         6
19  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      3.426831         6
20  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      0.913357         6
21  ApplyOfflineMap --in_data data/testCS-30-F1-O1...      0.358555         6
22  ApplyOfflineMap --in_data data/testCS-11-F2-O1...      0.171830         6
23  ApplyOfflineMap --in_data data/testCS-30-F1-O2...      0.797607         6
24  ApplyOfflineMap --in_data data/testCS-30-F1-O3...      1.484621         6
25  ApplyOfflineMap --in_data data/testCS-30-F1-O4...      2.175035         6
26  ApplyOfflineMap --in_data data/testCS-30-F2-O1...      0.276248         6
27  ApplyOfflineMap --in_data data/testCS-11-F1-O1...      0.073269         6
```

#### baseline/regtime.pkl
```
                                              command  average_time       time   per_diff
0   GenerateICOMesh --dual --res 60 --file meshes/...      0.528843   0.590863  11.727343
1   GenerateICOMesh --dual --res 30 --file meshes/...      0.185874   0.317657  70.898803
2   GenerateRLLMesh --lon 30 --lat 60 --file meshe...      0.107245   0.085638 -20.147216
3   GenerateCSMesh --file meshes/outCSMesh-30.g --...      0.085278   0.073856 -13.393983
4   GenerateCSMesh --file meshes/outCSMesh-11.g --...      0.074093   0.077373   4.426439
5   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.374052   0.348378  -6.863562
6   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.353323   0.379395   7.379042
7   GenerateTestData --mesh meshes/outCSMesh-11.g ...      0.125419   0.126770   1.077144
8   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.358127   0.353689  -1.239302
9   GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.353208   0.339167  -3.975133
10  GenerateTestData --mesh meshes/outCSMesh-30.g ...      0.364196   0.396182   8.782643
11  GenerateTestData --mesh meshes/outCSMesh-11.g ...      0.119653   0.111271  -7.004956
12  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...      0.400103   0.396702  -0.850087
13  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...      0.976406   0.994665   1.870063
14  GenerateOverlapMesh --a meshes/outCSMesh-30.g ...      4.320410   4.122126  -4.589466
15  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      2.049918   2.537829  23.801490
16  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      6.204872   6.895231  11.126073
17  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      4.293806   4.947679  15.228294
18  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      9.828979  10.430620   6.121086
19  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      0.788326   0.697766 -11.487678
20  GenerateOfflineMap --in_mesh meshes/outCSMesh-...      0.317520   0.276900 -12.793041
21  ApplyOfflineMap --in_data data/testCS-30-F1-O1...      0.250403   0.275291   9.939475
22  ApplyOfflineMap --in_data data/testCS-11-F2-O1...      0.106467   0.120101  12.805211
23  ApplyOfflineMap --in_data data/testCS-30-F1-O2...      0.574752   0.584282   1.658107
24  ApplyOfflineMap --in_data data/testCS-30-F1-O3...      1.084466   1.001624  -7.638974
25  ApplyOfflineMap --in_data data/testCS-30-F1-O4...      1.718680   1.746691   1.629851
26  ApplyOfflineMap --in_data data/testCS-30-F2-O1...      0.227098   0.194578 -14.319624
27  ApplyOfflineMap --in_data data/testCS-11-F1-O1...      0.067408   0.073388   8.871458
```