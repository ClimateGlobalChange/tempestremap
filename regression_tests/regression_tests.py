import subprocess
import os
import sys
import time
import errno
import pickle
import argparse
import unittest

from decimal import Decimal
import numpy as np
import pandas as pd

import multiprocessing as mp
from multiprocessing import Pool

# some global paths for storing files
meshes_path = "meshes/"
maps_path = "maps/"
data_path = "data/"
log_path = "logs/"
baseline_path = "baseline/"

# Global variable storing regression results
regression_dict = {'id': [], 'source': [], 'target': [], 'error': []}

timing_dict = {'command':[], 'average_time': [], 'num_runs': []}
time_pkl_file = baseline_path+'timing_data.pkl'

objects = []
generate_baseline = False

exitsTimeFile = True
#  Open timingfile
try:
    time_file = pickle.load(open(time_pkl_file, 'rb+'))

    with (open(time_pkl_file , "rb+")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    df_timingfile = objects[0]
    print("opened existing timing file: ", time_pkl_file)

except (OSError, EOFError, IOError) as e:

    exitsTimeFile = False
    print("opening new file")
    if not os.path.exists(os.path.dirname(time_pkl_file)):
        try:
            os.makedirs(os.path.dirname(time_pkl_file))
        except OSError as exec: 
            if exec.errno != errno.EEXIST:
                raise
    time_file = open(time_pkl_file, 'wb')


def populate_timing_data(cmd, t, loc):
    mycmd =""
    # expand the list to form command to execute
    for j in range(len(cmd)):
        # just to remove the bin path from the command that is run
        if(j==0):
            cmd[j] = cmd[j].split(loc)[1]
        mycmd+=cmd[j]+" "

    # timing file is created from the first time object read from timing file is empty
    if objects == []:
        timing_dict['command'].append(mycmd)
        timing_dict['average_time'].append(t)
        timing_dict['num_runs'].append(1)
    # use df_timingfile datafram read from previous timing
    elif generate_baseline == True:
        m = np.where(df_timingfile['command']==mycmd)
        val=(m[0]).item()

        # update the loaded file dataframe to put average time and number of runs for the command
        df_timingfile.iat[val, df_timingfile.columns.get_loc('num_runs')] = df_timingfile['num_runs'][val]+1
        df_timingfile.iat[val, df_timingfile.columns.get_loc('average_time')] = (df_timingfile['average_time'][val]+t)/2.0
    # no need to add timing to baselines, just compare and report differences  
    else:
        timing_dict['command'].append(mycmd)
        timing_dict['average_time'].append(t)
        timing_dict['num_runs'].append(1)

# a function to run a command and
# parse the output. 
def run_command(cmd):

    # using the Popen function to execute the command 
    mycmd =""
    # expand the list to form command to execute
    for j in range(len(cmd)):
        mycmd+=cmd[j]+" "

    cmd = mycmd
    # # use time command to measure the time take to run the task
    # cmd = "time "+ cmd
    if verbose:
        print("Running:", cmd,)
    temp = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
    timeStarted = time.time()         

    # TODO: dump each command line output to file(s)
    res = []
    success = True

    # use the communicate function to fetch the output
    output, errs = temp.communicate()
    timeDelta = time.time() - timeStarted                     # Get execution time.

    # parse pickle file for getting old counter and command
    counter = 0

    timing_dict['command'].append(cmd)
    timing_dict['average_time'].append(timeDelta)
    timing_dict['num_runs'].append(counter)


    # print("now timing: ", timing_dict)

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
      #  if we find the word EXCEPTION in the output report Failure
        if "EXCEPTION" in line:
           print("Error running ", cmd, "\n\n", line)
           success = False 
        res.append(line)



    return res, success, timeDelta

def to_float(inp):
    #return Decimal(inp)
    return np.float64(inp)

def extract_results(id, res):

    opts = []
    for i in range(1, len(res) - 1):
        line = res[i]

        #  look for the pattern below to find mass values
        a = line.find("....")
        if a != -1:
            opt = line.partition("....")[2].split()[2]
            # parse mass value for source and target both preceeded by .... eg. below value is 2.513274122871826e+01
            # ['Source', 'Mass:', '2.513274122871826e+01', 'Min', '1.0027367453e+00', 'Max', '2.9972632547e+00']
            opts.append(opt)
            # print(line.partition("....")[2].split()[0], " ", opt)

    #print('Found output: ', opts)
    srcV = to_float(opts[0])
    tgtV = to_float(opts[1])
    diff = abs(srcV - tgtV)
    percentage_diff = 100 * ( diff / srcV )

    regression_dict['id'].append(id)
    regression_dict['source'].append(srcV)
    regression_dict['target'].append(tgtV)
    regression_dict['error'].append(percentage_diff)

def generate_cs_mesh(res, filename):
    command = []
    command.append(bin_path+"/GenerateCSMesh")
    command.append("--file")
    command.append(filename)
    command.append("--res")
    command.append(res)

    return command

def generate_ico_mesh(res, dual, filename):
    command = []
    command.append(bin_path+"/GenerateICOMesh")
    if dual:
        command.append("--dual")
    command.append("--res")
    command.append(res)    
    command.append("--file")
    command.append(filename)

    return command

def generate_rll_mesh(lon, lat, filename):
    command = []
    command.append(bin_path+"/GenerateRLLMesh")
    command.append("--lon")
    command.append(lon)
    command.append("--lat")
    command.append(lat)
    command.append("--file")
    command.append(filename)

    return command

def generate_overlap_mesh(inpfname1, inpfname2, method, out):
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
def generate_offline_map(inpfname1, inpfname2, inpoverlapmesh, outputmap, orders, methods, correct_areas=False, monotone=False):
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
    
    if correct_areas:
        command.append("--correct_areas")
    
    if monotone:
        command.append("--mono")
    
    command.append("--out_map")
    command.append(maps_path+outputmap)    

    return command

def generate_test_data(inpfname, testname, out_test):
    command = []
    command.append(bin_path+"/GenerateTestData")
    command.append("--mesh")
    command.append(meshes_path+inpfname)
    command.append("--test")
    command.append(testname)
    command.append("--out")
    command.append(data_path+out_test)             

    return command

def apply_offline_map(mapfile, inputdatafile, variablename, outfile):
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

def generate_mesh(meshstr):
    if "cs" in meshstr:
        filename = "outCSMesh"
        res = meshstr.partition("cs-")[2]
        filename+="-"+res+".g"

        # call function to generate cs mesh command
        command = generate_cs_mesh(res, meshes_path+filename)
        #
    elif "icod" in meshstr:
        filename = "outICODMesh"
        res = meshstr.partition("icod-")[2]
        filename+="-"+res+".g"

        # call function to generate polygonal mesh command        
        command = generate_ico_mesh(res, True, meshes_path+filename)
        #
    elif "ico" in meshstr:
        filename = "outICOMesh"
        res = meshstr.partition("ico-")[2]
        filename+="-"+res+".g"

        # call function to generate icosahedral mesh command        
        command = generate_ico_mesh(res, False, meshes_path+filename)
        #
    elif "rll" in meshstr:
        filename = "outRLLMesh"
        tres = meshstr.partition("rll-")[2]
        lon = tres.partition("-")[0]
        lat = tres.partition("-")[2]
        filename+="-"+lon+"-"+lat+".g"

        # call function to generate regular lat-lon mesh command
        command = generate_rll_mesh(lon, lat, meshes_path+filename)
        #

    return command, filename

def check_error(success):
    if False in success:
       print("exiting..")
       print("Find usage by running: \npython regression_tests.py -h\n")

       exit()   

if __name__ == '__main__':

    assert sys.version_info >= (3,0) 

    # default arguments
    verbose=False
    procs = 2

    parser=argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-g", "--generate_baseline", help="Generate new baseline for comparison (regression) later", action = "store_true")
    parser.add_argument("-p", "--path", type=str, help="run parallel tests", default="../bin")
    parser.add_argument("-n", "--procs", type=int, help="number of processors used for running regression tests",default=2)

    # ex arguments:
    # generate baselines:                   python regression_tests.py -v -g -n 3 -p ../bin/ 
    # regression compare against baselines: python regression_tests.py -v -n 3 -p ../bin/ 
    args = parser.parse_args()
    if args.verbose:
        print("verbose mode enabled")
        verbose=True
    if args.procs:
        print("tasks specified=", args.procs)
        if args.procs > 0: procs=args.procs
    if args.generate_baseline:
        print("Generating baseline")
        generate_baseline=True
    if args.path:
        print("executable path specified =", args.path)
        bin_path=args.path


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

    # this is used to remove the path from the command later for storing timing results
    splitLoc = bin_path + "/"

    success = []
    results = []
    count = 0
    #  print(mesh_cmds)

    # create mesh data maps and test directories before running if not already created
    if not os.path.isdir("meshes"):
        os.mkdir("meshes")
    if not os.path.isdir("data"):
        os.mkdir("data")    
    if not os.path.isdir("test"):
        os.mkdir("test")
    if not os.path.isdir("maps"):
        os.mkdir("maps")

    # Each for loop below runs the commands to create results
    for result, ss, timeDelta in pool.map(run_command, mesh_cmds):

        populate_timing_data(mesh_cmds[count], timeDelta, splitLoc)
        count=count+1

        results.append(result)
        success.append(ss)
    check_error(success)
    #  
    success = []
    results = []
    count =0
    #  print(generate_test_cmds) 
    for result, ss, timeDelta in pool.map(run_command, generate_test_cmds):

       populate_timing_data(generate_test_cmds[count], timeDelta, splitLoc)
       count=count+1
       results.append(result)
       success.append(ss)    
    check_error(success)
    #
    success = []
    results = [] 
    count = 0
    #  print(overlap_test_cmds)
    for result, ss, timeDelta in pool.map(run_command, overlap_test_cmds):

       populate_timing_data(overlap_test_cmds[count], timeDelta, splitLoc)
       count=count+1
       results.append(result)
       success.append(ss)    
    check_error(success)
    # 
    success = []
    results = []    
    count = 0
    #  print(g_offmap_cmds)
    for result, ss, timeDelta in pool.map(run_command, g_offmap_cmds):

       populate_timing_data(g_offmap_cmds[count], timeDelta, splitLoc)
       count=count+1
       results.append(result)
       success.append(ss)    
    check_error(success)
    # 
    success = []
    results = []   
    count =0 
    #  print(a_offmap_cmds)
    for result, ss, timeDelta in pool.map(run_command, a_offmap_cmds):

       populate_timing_data(a_offmap_cmds[count], timeDelta, splitLoc)
       count=count+1
       results.append(result)
       success.append(ss)       
    check_error(success)
    #  print(results)
    pool.close()
    pool.join()

    pd.options.display.float_format = '{:,2.16f}'.format

    # create baseline results
    if generate_baseline:
        count = 0
        for i in range(len(tm.values)):
            id = tm.loc[i].at["id"]
            extract_results(id, results[count])
            count+=1

        df = pd.DataFrame(regression_dict)
        df['source'] = df['source'].astype(np.float64)
        df['target'] = df['target'].astype(np.float64)
        df['error'] = df['error'].astype(np.float64)

        # saving the dataframe
        base_pkl_file = baseline_path+'baseline_data.pkl'
        pkl_file = open(base_pkl_file, 'wb')
        pickle.dump(df, pkl_file)
        print("\n Saved baseline/baseline_data.pkl file")

        # write current results to a file
        #  check if objects populated from existing timing file is available
        if objects == []:
            df_t = pd.DataFrame(timing_dict)
            pickle.dump(df_t, time_file)
            print("\n Saved baseline/timing_data.pkl file")
        else:
            # open the file again with write mode
            with open(time_pkl_file, 'wb') as time_file:
                pickle.dump(df_timingfile, time_file)  
                print("\n Updated baseline/timing_data.pkl file")
    # compare against existing baseline
    else:
        # read baseline results from baseline_data.csv file
        base_pkl_file = baseline_path+'baseline_data.pkl'
        bpfile = open(base_pkl_file, 'rb')

        df = pickle.load(bpfile)
        df = df.set_index(['id'])
        df['source'] = df['source'].astype(np.float64)
        df['target'] = df['target'].astype(np.float64)
        df['error'] = df['error'].astype(np.float64)

        count = 0
        print('Difference against baselines\n                      id, source-diff target-diff error-diff')

        assert df.shape[0] >= len(tm.values), "baseline has fewer number of items than current run\n generate baselines first.."

        # extract current results
        for i in range(len(tm.values)):
            id = tm.loc[i].at["id"]
            extract_results(id, results[count])
            df_current = pd.DataFrame(regression_dict)
            df_current = df_current.set_index(['id'])
            df_current['source'] = df_current['source'].astype(np.float64)
            df_current['target'] = df_current['target'].astype(np.float64)
            df_current['error'] = df_current['error'].astype(np.float64)

            print('Baseline comparison: ', id, df['source'][id]-df_current['source'][id], 
                                              df['target'][id]-df_current['target'][id], 
                                              df['error'][id]-df_current['error'][id])
            count += 1

        # write current results to a file
        rpt_pkl_file = baseline_path+'baseline_data_repeat.pkl'
        pkl_file = open(rpt_pkl_file, 'wb')
        pickle.dump(df_current, pkl_file)
        print("\n Saved baseline/baseline_data_repeat.pkl file")

        # write a new timing pickle file that compares against the previous run
        # cmd avg/time_baseline  current_average_time  percentage_difference
        # cmd, num_baseline_runs, baseline_time, time, per_diff
        comptime_dict = {'command':[], 'num_runs':[], 'average_time': [], 'time': [], 'per_diff':[]}
        comptime_dict = {'command':[], 'average_time': [], 'time': [], 'per_diff':[]}

        # timing_dict has current timing  
        # df_timingfile has results from old file.
        # print(timing_dict['average_time'])
        dlength=len(timing_dict['average_time'])
        # print(dlength)
        # print(df_timingfile.shape[0], df_timingfile.shape[1])
        # print("old_results")
        print(df_timingfile['average_time'][0])
        for j in range(dlength):
            m = np.where(df_timingfile['command']==timing_dict['command'][j])
            # assuming there are no repetitions and there is one unique cmd
            val=np.asscalar(m[0])
            # print(val)
            comptime_dict['command'].append(timing_dict['command'][j])
            # comptime_dict['num_runs'].append(df_timingfile['num_runs'][val])
            comptime_dict['average_time'].append(df_timingfile['average_time'][val])
            comptime_dict['time'].append(timing_dict['average_time'][j])
            
            per_diff = (timing_dict['average_time'][j] - df_timingfile['average_time'][val])*100/df_timingfile['average_time'][val]

            comptime_dict['per_diff'].append(per_diff)
        
        
        dfcomp = pd.DataFrame(comptime_dict)
        dfcomp['command'] = dfcomp['command'].astype(np.str)
        # dfcomp['num_runs'] = dfcomp['num_runs'].astype(np.int)
        dfcomp['average_time'] = dfcomp['average_time'].astype(np.float64)
        dfcomp['time'] = dfcomp['time'].astype(np.float64)
        dfcomp['per_diff'] = dfcomp['per_diff'].astype(np.float64)


        # saving the dataframe
        regtime_pkl_file = baseline_path+'regtime.pkl'
        regt_file = open(regtime_pkl_file, 'wb')
        pickle.dump(dfcomp, regt_file)
        print("\n Saved baseline/regtime.pkl file")

        # print(comptime_dict)