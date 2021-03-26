#coding=utf8

########################################################################

def get_nodes_from_node_set(mesh_file_name, node_set_name):

    mesh_file = open(mesh_file_name, "r")

    nodes_all = {}
    nodes_sel = []

    context = ""
    for line in mesh_file:
        if (line[-1:] == "\n"): line = line[:-1]
        #if (verbose): print "line =", line

        if line.startswith("**"): continue

        if (context == "reading nodes"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                nodes_all[int(splitted_line[0])] = [float(coord) for coord in splitted_line[1:4]]

        if (context == "reading node set"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                nodes_sel += [int(node_number) for node_number in splitted_line if node_number != '']

        if (context == "reading nodes and node set"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                nodes_all[int(splitted_line[0])] = [float(coord) for coord in splitted_line[1:4]]
                nodes_sel += [int(splitted_line[0])]

        if line.startswith("*NODE"):
            if ("NSET="+node_set_name in line):
                context = "reading nodes and node set"
            else:
                context = "reading nodes"
        if line.startswith("*NSET") and ("NSET="+node_set_name in line):
            context = "reading node set"

    mesh_file.close()

    return {node_num:node_pos for node_num,node_pos in nodes_all.iteritems() if node_num in nodes_sel}
