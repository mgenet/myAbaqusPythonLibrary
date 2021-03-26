#coding=utf8

########################################################################

def get_elems_from_elem_set(mesh_file_name, elem_set_name):

    mesh_file = open(mesh_file_name, "r")

    elems_all = {}
    elems_sel = []

    context = ""
    for line in mesh_file:
        if (line[-1:] == "\n"): line = line[:-1]
        #if (verbose): print "line =", line

        if line.startswith("**"): continue

        if (context == "reading elems"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                elems_all[int(splitted_line[0])] = [int(node_number) for node_number in splitted_line[1:]]

        if (context == "reading elem set"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                elems_sel += [int(elem_number) for elem_number in splitted_line if elem_number != '']

        if (context == "reading elems and elem set"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                elems_all[int(splitted_line[0])] = [int(node_number) for node_number in splitted_line[1:]]
                elems_sel += [int(splitted_line[0])]

        if line.startswith("*ELEMENT"):
            if ("ELSET="+elem_set_name in line):
                context = "reading elems and elem set"
            else:
                context = "reading elems"
        if line.startswith("*ELSET") and ("ELSET="+elem_set_name in line):
            context = "reading elem set"

    mesh_file.close()

    return {elem_num:elem_nodes for elem_num,elem_nodes in elems_all.items() if elem_num in elems_sel}
