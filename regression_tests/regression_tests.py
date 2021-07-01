import subprocess
import os
import pandas as pd
import sys
import errno

import multiprocessing as mp
from multiprocessing import Pool

bin_path = "../bin/"

# a function to run a command and
# parse the output.
def run_command(cmd):

    # using the Popen function to execute the command 
    mycmd =""
    # expand the list to form command to execute
    for j in range(len(cmd)):
        mycmd+=cmd[j]+" "

    cmd = mycmd     
    print("cmd running:", cmd,)
    temp = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)

    # TODO: dump each command line output to file(s)
    res = []
    success = True

    # use the communicate function to fetch the output
    output, errs = temp.communicate()
   
   # if the output is empty - executables not found then return False   
    if output == b'':
       success = False

    # splitting the output to parse line by line
    output = str(output)
    output = output.split("\n")

    output = output[0].split('\\n')

    # a variable to store the output
    res = []

    # iterate line by line
    for line in output:
      #   print (line)
      #  if we find the word EXCEPTION in the output report Failure
        if "EXCEPTION" in line:
           print("Error running ", cmd, "\n\n", line)
           success = False 
        res.append(line)

    return res, success

def create_baseline_results(id, res):

    opts = []

    filename = "./baseline/regression_results.txt"
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    with open(filename, "a") as f:
         for i in range(1, len(res) - 1):
               line = res[i]
               # if "--in_data <string>" in line or "--map <string> " in line:
               #    print(line)
               #    # f.write(line)
               #    # f.write("\n")

               #  look for the pattern below to find mass values
               a = line.find("....")
               if a != -1:
                  opt = line.partition("....")[2].split()[2]
                  # parse mass value for source and target both preceeded by .... eg. below value is 2.513274122871826e+01
                  # ['Source', 'Mass:', '2.513274122871826e+01', 'Min', '1.0027367453e+00', 'Max', '2.9972632547e+00']
                  opts.append(opt)
                  # print(line.partition("....")[2].split()[0], " ", opt)


         diff = float(opts[0]) - float(opts[1])
         percentage_diff = ( diff * 100 )/float(opts[0])
         # f.write("\n")
         text = "CaseID: " + str(id) + ": src val, target val, err " 
         f.write("\n")
         f.write(text)
         f.write("\n")
         vals =  str(id) + " " + str(opts[0]) + " " + str(opts[1]) + " " + str(percentage_diff)
         f.write(vals)
         f.write("\n")

         print(vals)   

def generate_cs_mesh(res, filename):
    command = []
    command.append(bin_path+"GenerateCSMesh")
    command.append("--file")
    command.append(filename)
    command.append("--res")
    command.append(res)

    return command

def generate_ico_mesh(res, filename, dual=False):
    command = []
    command.append(bin_path+"GenerateICOMesh")
    command.append("--dual")
    # command.append(dual)
    command.append("--res")
    command.append(res)    
    command.append("--file")
    command.append(filename)

    return command

def generate_rll_mesh(lon, lat, filename):
    command = []
    command.append(bin_path+"GenerateRLLMesh")
    command.append("--lon")
    command.append(lon)
    command.append("--lat")
    command.append(lat)
    command.append("--file")
    command.append(filename)

    return command

def generate_overlap_mesh(inpfname1, inpfname2, method, out):
    command = []
    command.append(bin_path+"GenerateOverlapMesh")
    command.append("--a")
    command.append("meshes/"+inpfname1)
    command.append("--b")
    command.append("meshes/"+inpfname2)
    command.append("--out")
    command.append("meshes/"+out)    
    command.append("--method")
    command.append(method)

    return command

# order: list with 2 entries that is strtictly above 0: FV (1:4), SE==4
# methods: list with 2 entries specifying one of fv, cgll, dgll
def generate_offline_map(inpfname1, inpfname2, inpoverlapmesh, outputmap, orders, methods, correct_areas=False, monotone=False):
    command = []
    command.append(bin_path+"GenerateOfflineMap")
    command.append("--in_mesh")
    command.append("meshes/"+inpfname1)
    command.append("--out_mesh")
    command.append("meshes/"+inpfname2)
    command.append("--ov_mesh")
    command.append("meshes/"+inpoverlapmesh)
    command.append("--in_np")
    command.append(orders[0])
    command.append("--out_np")
    command.append(orders[1])
    command.append("--in_type")
    command.append(methods[0])
    command.append("--out_type")
    command.append(methods[1])
    command.append("--correct_areas")
    # command.append(correct_areas)
    command.append("--mono")
    # command.append(monotone)
    command.append("--out_map")
    command.append("maps/"+outputmap)    

    return command

def generate_test_data(inpfname, testname, out_test):
    command = []
    command.append(bin_path+"GenerateTestData")
    command.append("--mesh")
    command.append("meshes/"+inpfname)
    command.append("--test")
    command.append(testname)
    command.append("--out")
    command.append("data/"+out_test)             

    return command

def apply_offline_map(mapfile, inputdatafile, variablename, outfile):
    command = []
    command.append(bin_path+"ApplyOfflineMap")
    command.append("--in_data")
    command.append("data/"+inputdatafile)
    command.append("--map")
    command.append("maps/"+mapfile)
    command.append("--var")
    command.append(variablename)    
    command.append("--out_data")
    command.append("data/"+outfile)    
        
    return command

def generate_mesh(meshstr):
    if "cs" in meshstr:
        filename = "outCSMesh"
        res = meshstr.partition("cs-")[2]
        filename+="-"+res+".g"

        # call function to generate cs mesh command
        command = generate_cs_mesh(res, "meshes/"+filename)
        #
    elif "icod" in meshstr:
        filename = "outICOMesh"
        res = meshstr.partition("icod-")[2]
        filename+="-"+res+".g"

        # call function to generate cs mesh command        
        command = generate_ico_mesh(res, "meshes/"+filename)
        #
    elif "rll" in meshstr:
        filename = "outRLLMesh"
        tres = meshstr.partition("rll-")[2]
        lon = tres.partition("-")[0]
        lat = tres.partition("-")[2]
        filename+="-"+lon+"-"+lat+".g"

        # call function to generate cs mesh command
        command = generate_rll_mesh(lon, lat, "meshes/"+filename)
        #

    return command, filename

def check_error(success):
    if False in success:
       print("exiting..")
       print("Usage: \npython regression_tests.py <location of executables>  <number of procs>\n    Default <location of executables> is ../bin/ and <number of procs> is 2\n")

       exit()   

if __name__ == '__main__':
    procs = 2
    print("Usage: \npython regression_tests.py <location of executables>  <number of procs>\n    Default <location of executables> is ../bin/ and <number of procs> is 2\n")
    if len(sys.argv) >= 2:
       bin_path = sys.argv[1] + "/"

    if len(sys.argv) == 3:
        procs= int (sys.argv[2])

    # Run a pipeline
    # read inputs
    command = []
    args=[]
    # space seperated table with keywords
    tm = pd.read_table('./test_matrix.ini', delim_whitespace=True, comment='#')

    # figure out order of commands to call
    # For each pipeline:
    # this loop gets the commands to call and the arguments, size of both commands and arguments is same

    # collect all mesh generation commands
    mesh_cmds = []
    generate_test_cmds = []
    overlap_test_cmds = []
    g_offmap_cmds = []
    a_offmap_cmds = []

    for i in range(len(tm.values)):

        print("Running test case ",  tm.loc[i].at["id"])
# Generate Meshes
        srcmesh_str =  tm.loc[i].at["srcmesh"]
        tgtmesh_str = tm.loc[i].at["tgtmesh"]

        cmd, in_file = generate_mesh(srcmesh_str)
        mesh_cmds.append(cmd)
        
        cmd, out_file = generate_mesh(tgtmesh_str)
        mesh_cmds.append(cmd)       

# Generate Overlap Mesh
        method = "exact"

        # call overlap mesh
        # name overlap mesh as in+outfile.g
        overlapmesh = in_file[:len(in_file)-2]+"-"+out_file

        cmd = generate_overlap_mesh(in_file, out_file, method, overlapmesh)
        overlap_test_cmds.append(cmd)

# Generate Offline Map
        order = str(tm.loc[i].at["order"])
        ofmap_out="mapNE-"
        ofmap_out+=srcmesh_str.upper()+"-"+tgtmesh_str.upper()
        ofmap_out+="-O"+order+".nc"

        correct_areas = tm.loc[i].at["correctareas"]
        monotone = tm.loc[i].at["monotone"]

        # specify two entries 1st for input and 2nd for output
        orders = []
        methods = []

        orders.append(order)
        orders.append(order)

        methods.append("fv")
        methods.append("fv")

        inpoverlapmesh = overlapmesh 

        # call generate offline map
        cmd = generate_offline_map(in_file, out_file, inpoverlapmesh, ofmap_out, orders, methods, correct_areas, monotone)
        g_offmap_cmds.append(cmd)

# Generate Test Data
        test = str(tm.loc[i].at["test"])
        out_td_src="test"
        #srcmesh
        out_td_src+=srcmesh_str.upper()
        out_td_src+="-F"+test+"-O"+order+".nc"

        cmd = generate_test_data(in_file, test, out_td_src)
        generate_test_cmds.append(cmd)

      #   out_td_tgt="test"
      #   #tgtmesh
      #   out_td_tgt+=tgtmesh_str.upper()
      #   out_td_tgt+="-F"+test+"-O"+order+".nc"
      #   cmd = generate_test_data(out_file, test, out_td_tgt)
      #   generate_test_cmds.append(cmd)

# Apply Offline Map
        variablename = "Psi"
        out_test="OutTest"
        out_test+=srcmesh_str.upper()
        out_test+="-F"+test+"-O"+order+".nc"

        cmd = apply_offline_map(ofmap_out, out_td_src, variablename, out_test)          
        a_offmap_cmds.append(cmd)  

    

# Running multiple processes in parallel and in sequence
    pool = Pool(processes=procs)
   #  
   # get unique commands to run, convert list to set and back
    mesh_cmds = [list(x) for x in set(tuple(x) for x in mesh_cmds)] 
    generate_test_cmds = [list(x) for x in set(tuple(x) for x in generate_test_cmds)] 
    overlap_test_cmds = [list(x) for x in set(tuple(x) for x in overlap_test_cmds)] 
    g_offmap_cmds = [list(x) for x in set(tuple(x) for x in g_offmap_cmds)] 
    a_offmap_cmds = [list(x) for x in set(tuple(x) for x in a_offmap_cmds)]
   # 

    success = []
    results = []
    count = 0
   #  print(mesh_cmds)
    for result, ss in pool.map(run_command, mesh_cmds):
       results.append(result)
       success.append(ss)
    check_error(success)
   #  
    success = []
    results = []
   #  print(generate_test_cmds) 
    for result, ss in pool.map(run_command, generate_test_cmds):
       results.append(result)
       success.append(ss)    
    check_error(success)
   #
    success = []
    results = [] 
   #  print(overlap_test_cmds)
    for result, ss in pool.map(run_command, overlap_test_cmds):
       results.append(result)
       success.append(ss)    
    check_error(success)
   # 
    success = []
    results = []    
   #  print(g_offmap_cmds)
    for result, ss in pool.map(run_command, g_offmap_cmds):
       results.append(result)
       success.append(ss)    
    check_error(success)
   # 
    success = []
    results = []    
   #  print(a_offmap_cmds)
    for result, ss in pool.map(run_command, a_offmap_cmds):
       results.append(result)
       success.append(ss)       
    check_error(success)
   #  print(results)
    pool.close()
    pool.join()

    # create baseline results
    count = 0
    for i in range(len(tm.values)):
      id = tm.loc[i].at["id"]
      create_baseline_results(id, results[count])
      count+=1

    # TODO: compare new against baselines



