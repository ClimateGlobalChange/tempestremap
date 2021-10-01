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
                                              command  num_runs  average_time       time    per_diff
0   GenerateCSMesh --file meshes/outCSMesh-30.g --...         5      0.087082   0.223496  156.651384
1   GenerateRLLMesh --lon 30 --lat 60 --file meshe...         5      0.086561   0.212475  145.463647
2   GenerateICOMesh --dual --res 30 --file meshes/...         5      0.277692   0.261011   -6.007078
3   GenerateICOMesh --dual --res 60 --file meshes/...         5      0.717480   0.613342  -14.514360
4   GenerateCSMesh --file meshes/outCSMesh-11.g --...         5      0.059156   0.048707  -17.663702
5   GenerateTestData --mesh meshes/outCSMesh-11.g ...         5      0.110063   0.155130   40.947451
6   GenerateTestData --mesh meshes/outCSMesh-30.g ...         5      0.519196   0.398670  -23.213861
7   GenerateTestData --mesh meshes/outCSMesh-30.g ...         5      0.463584   0.470161    1.418859
8   GenerateTestData --mesh meshes/outCSMesh-30.g ...         5      0.432348   0.460484    6.507617
9   GenerateTestData --mesh meshes/outCSMesh-30.g ...         5      0.451841   0.530935   17.504841
10  GenerateTestData --mesh meshes/outCSMesh-11.g ...         5      0.121148   0.115940   -4.298604
11  GenerateTestData --mesh meshes/outCSMesh-30.g ...         5      0.485198   0.330369  -31.910593
12  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...         5      0.562449   0.546039   -2.917613
13  GenerateOverlapMesh --a meshes/outCSMesh-11.g ...         5      1.948958   1.901371   -2.441649
14  GenerateOverlapMesh --a meshes/outCSMesh-30.g ...         5      6.236995   5.923357   -5.028666
15  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5     22.833240  25.009068    9.529213
16  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5      1.096764   1.135980    3.575616
17  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5      0.388046   0.367089   -5.400683
18  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5     10.007985  10.073919    0.658807
19  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5     14.706096  15.428331    4.911131
20  GenerateOfflineMap --in_mesh meshes/outCSMesh-...         5      4.613005   3.326606  -27.886356
21  ApplyOfflineMap --in_data data/testCS-30-F1-O1...         5      0.263598   0.492674   86.903816
22  ApplyOfflineMap --in_data data/testCS-11-F2-O1...         5      0.121074   0.210612   73.953442
23  ApplyOfflineMap --in_data data/testCS-30-F1-O2...         5      0.947658   1.715561   81.031617
24  ApplyOfflineMap --in_data data/testCS-30-F1-O3...         5      1.956097   2.706602   38.367462
25  ApplyOfflineMap --in_data data/testCS-30-F1-O4...         5      2.604080   2.677558    2.821637
26  ApplyOfflineMap --in_data data/testCS-30-F2-O1...         5      0.317955   0.310593   -2.315633
27  ApplyOfflineMap --in_data data/testCS-11-F1-O1...         5      0.077363   0.329290  325.643913
```