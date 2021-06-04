import subprocess
import os
import pandas as pd

# a function to run a command and
# parse the output.
def run_command(cmd):

    # using the Popen function to execute the
    # command and store the result in temp.
    # it returns a tuple that contains the
    # data and the error if any.
    print("cmd running:", cmd,)
    temp = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)

    # we use the communicate function
    # to fetch the output
    output = str(temp.communicate())

    # splitting the output so that
    # we can parse them line by line
    output = output.split("\n")

    output = output[0].split('\\')

    # a variable to store the output
    res = []

    # iterate through the output
    # line by line
    for line in output:
        res.append(line)

    return res

def parse_help(res):
    opts = []
    for i in range(1, len(res) - 1):
      line = res[i]
      #print("line", i, line)
      a = line.find("--")
      if a != -1:
        opt = line.partition("--")[2].split()[0]
        opts.append(opt)
        print("line is", opt)

    return opts

def create_results(res):
    opts = []
    for i in range(1, len(res) - 1):
      line = res[i]
      #print("line", i, line)
      a = line.find("....")
      if a != -1:
        opt = line.partition("....")[2].split()[2]
# ['Source', 'Mass:', '2.513274122871826e+01', 'Min', '1.0027367453e+00', 'Max', '2.9972632547e+00']
        opts.append(opt)
        print(line.partition("....")[2].split()[0], " ", opt)

    return opts

def build_mesh_args(value, type):
    res = value
    args=[]
    args = " --res "+res
    if type == "cs":
        filename = "outCSMesh"+res+".g"
        args+= " --file "+ filename
    elif type == "icod":
        filename = "outICOMesh"+res+".g"
        args+= " --file "+ filename

    return args, filename

if __name__ == '__main__':

    # These option are availble for checking and extending this regression testing commandline options
    res = run_command('./GenerateCSMesh')
    gCSMesh_tokens = parse_help(res)
    res = run_command('./GenerateICOMesh')
    gICOMesh_tokens = parse_help(res)
    res = run_command('./GenerateOverlapMesh')
    gOverlapMesh_tokens = parse_help(res)
    res = run_command('./GenerateOfflineMap')
    gOfflineMap_tokens = parse_help(res)
    res = run_command('./GenerateTestData')
    gTestData_tokens = parse_help(res)
    res = run_command('./ApplyOfflineMap')
    gApplyOfflineMap_tokens = parse_help(res)

    # Run a pipeline
    # read inputs
    command = []
    args=[]
    # delim_whitespace=True is needed is text in csv file is seperated by just spaces
    tm = pd.read_csv('./test_matrix.csv', delim_whitespace=True)
    # figure out order of commands to call
    for i in range(len(tm.values)):
# Command 1
        if "cs" in tm.loc[i].at["srcmesh"]:
            # we have a CSMesh
            command.append("GenerateCSMesh")
    # Finalize parameters for the command
            srcvalue = tm.loc[i].at["srcmesh"].partition("cs")[2]
            arg, in_file = build_mesh_args(srcvalue, "cs")
            args.append(arg)
            flag_a = "cs"
        elif "icod" in tm.loc[i].at["srcmesh"]:
            command.append("GenerateICOMesh")
            srcvalue = tm.loc[i].at["srcmesh"].partition("icod")[2]
            arg, in_file = build_mesh_args(srcvalue, "icod")
            arg+= " --dual "
            args.append(arg)
            flag_a = "icod"

# Command 2
        if "cs" in tm.loc[i].at["tgtmesh"]:
            # we have a CSMesh
            command.append("GenerateCSMesh")
            tgtvalue = tm.loc[i].at["tgtmesh"].partition("cs")[2]
            arg, out_file = build_mesh_args(tgtvalue, "cs")
            args.append(arg)
            flag_b = "cs"

        elif "icod" in tm.loc[i].at["tgtmesh"]:
            # we have a ICOMesh
            command.append("GenerateICOMesh")
            tgtvalue = tm.loc[i].at["tgtmesh"].partition("icod")[2]
            arg, out_file = build_mesh_args(tgtvalue, "icod")
            arg+= " --dual "
            args.append(arg)
            flag_b = "icod"

# Command 3
        command.append("GenerateOverlapMesh")
        arg = " --a "+ in_file
        arg+= " --b " + out_file
        arg+= " --method exact"

        args.append(arg)

# Command 4
        command.append("GenerateOfflineMap")
        arg = " --ov_mesh overlap.g"
        order = str(tm.loc[i].at["order"])
        arg+= " --in_np "+ order+ " --out_np "+ order
        arg+= " --out_mesh " + out_file+ " --in_mesh "+ in_file
        ofmap_out="mapNE-"
        #srcmesh
        if flag_a == "cs":
            ofmap_out+= "CS"+srcvalue
        elif flag_a == "icod":
            ofmap_out+= "ICOD"+srcvalue
        #tgtmesh
        if flag_b == "cs":
            ofmap_out+= "-CS"+tgtvalue
        elif flag_b == "icod":
            ofmap_out+= "-ICOD"+tgtvalue
        ofmap_out+="-O"+order+".nc"
        arg+= " --out_map " + ofmap_out

        args.append(arg)

# Command 5
        command.append("GenerateTestData")
        arg = " --mesh "+ in_file
        test = str(tm.loc[i].at["test"])
        arg+= " --test "+ test
        out_td="test"
        #srcmesh
        if flag_a == "cs":
            out_td+= "CS"+srcvalue
        elif flag_a == "icod":
            out_td+= "ICOD"+srcvalue

        out_td+="-F"+test+"-O"+order+".nc"
        arg+= " --out " + out_td

        args.append(arg)

# Command 6
        command.append("ApplyOfflineMap")
        arg = " --in_data "+ out_td
        arg+= " --map "+ ofmap_out
        arg+= " --var Psi"
        out_test="OutTest"
        #srcmesh
        if flag_a == "cs":
            out_test+= "CS"+srcvalue
        elif flag_a == "icod":
            out_test+= "ICOD"+srcvalue

        out_test+="-F"+test+"-O"+order+".nc"
        arg+= " --out_data " + out_test

        args.append(arg)




        tm.loc[i].at["correctareas"]
        tm.loc[i].at["monotone"]

    print("args =",args)
    print("commands:", command)
    # loop thru all commands and run
    # check if the same command has already been run
    path = "../build/"
    num_commands = len(command)
    res = []
    for i in range(num_commands):
        cmd=path+command[i]
        cmd_all = cmd+args[i]
        res.append(run_command(cmd_all))
        if command[i] == "ApplyOfflineMap":
            results = create_results(res[i])
            diff = float(results[0]) - float(results[1])
            percentage_diff = ( diff * 100 )/float(results[0])
            print("Percentage diff between source and target is ", percentage_diff)

    #  # make args for call
    #  # launch the commands
    # store output and compareds = []

