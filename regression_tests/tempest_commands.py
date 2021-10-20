def generate_cs_mesh(res, filename, bin_path):
    command = []
    command.append(bin_path+"/GenerateCSMesh")
    command.append("--file")
    command.append(filename)
    command.append("--res")
    command.append(res)

    return command

def generate_ico_mesh(res, dual, filename, bin_path):
    command = []
    command.append(bin_path+"/GenerateICOMesh")
    if dual:
        command.append("--dual")
    command.append("--res")
    command.append(res)    
    command.append("--file")
    command.append(filename)

    return command

def generate_rll_mesh(lon, lat, filename, bin_path):
    command = []
    command.append(bin_path+"/GenerateRLLMesh")
    command.append("--lon")
    command.append(lon)
    command.append("--lat")
    command.append(lat)
    command.append("--file")
    command.append(filename)

    return command

def generate_overlap_mesh(inpfname1, inpfname2, method, out, bin_path, meshes_path):
    command = []
    command.append(bin_path+"/GenerateOverlapMesh")
    command.append("--a")
    command.append(meshes_path+inpfname1)
    command.append("--b")
    command.append(meshes_path+inpfname2)
    command.append("--out")
    command.append(meshes_path+out)    
    command.append("--method")
    command.append(method)

    return command

# order: list with 2 entries that is strtictly above 0: FV (1:4), SE==4
# methods: list with 2 entries specifying one of fv, cgll, dgll
def generate_offline_map(inpfname1, inpfname2, inpoverlapmesh, outputmap, orders, methods, bin_path, meshes_path, maps_path, nocorrect_areas=False, monotone=False):
    command = []
    command.append(bin_path+"/GenerateOfflineMap")
    
    command.append("--in_mesh")
    command.append(meshes_path+inpfname1)
    command.append("--out_mesh")
    command.append(meshes_path+inpfname2)
    command.append("--ov_mesh")
    command.append(meshes_path+inpoverlapmesh)
    
    command.append("--in_np")
    command.append(orders[0])
    command.append("--out_np")
    command.append(orders[1])
   
    command.append("--in_type")
    command.append(methods[0])
    command.append("--out_type")
    command.append(methods[1])
    
    if nocorrect_areas:
        command.append("--nocorrectareas")
    
    if monotone:
        command.append("--mono")
    
    command.append("--out_map")
    command.append(maps_path+outputmap)    

    return command

def generate_test_data(inpfname, testname, out_test, bin_path, data_path, meshes_path):
    command = []
    command.append(bin_path+"/GenerateTestData")
    command.append("--mesh")
    command.append(meshes_path+inpfname)
    command.append("--test")
    command.append(testname)
    command.append("--out")
    command.append(data_path+out_test)             

    return command

def apply_offline_map(mapfile, inputdatafile, variablename, outfile, bin_path, data_path, maps_path):
    command = []
    command.append(bin_path+"/ApplyOfflineMap")
    
    command.append("--in_data")
    command.append(data_path+inputdatafile)
    
    command.append("--map")
    command.append(maps_path+mapfile)
    
    command.append("--var")
    command.append(variablename)    
    
    command.append("--out_data")
    command.append(data_path+outfile)    
        
    return command