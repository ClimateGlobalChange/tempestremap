***
Overview:
- Paths used for regression are relative to the directory <yourpath>/tempestremap/regression_tests/
    - Python calls must be made from this directory
    - There are 5 major commands in the regression workflow
      - Generate mesh (called twice)
        - GenerateICOMesh
        - GenerateCSMesh
        - GenerateRLLMesh
      - GenerateTestData
      - GenerateOfflineMap
      - GenerateOverlapMesh
      - ApplyOfflineMap
    - The test_matrix.ini file  specifies the command line options for each of the 6 commands
    - % difference b/w source and target is calculated and reported

More options and commands can be added

takes one argument that specifies the path to executables


*** 
Directory Structure
- regression_tests
    - baseline: pickle files for baseline runs, baseline timing and regression time.
    - data: test data files
    - maps: map files
    - meshes: mesh files

