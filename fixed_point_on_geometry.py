#coding=utf8

### ---------------------------------------------------------------- ###

import os
import sys
import numpy
import shutil

from wait_for_job_to_complete import *

### ---------------------------------------------------------------- ###

def fixed_point_on_geometry(sedlines,
                            reference_nodes_file_basename="Nodes",
                            initial_nodes_file_basename="Nodes",
                            input_file_basename="Job",
                            cpus=None,
                            gpus=None,
                            user=None,
                            tol=1.1e-3,
                            nb_max_iterations=100,
                            verbose=True):

    if (verbose): print "Read reference coordinates"
    X_ref = numpy.array([[float(coordinate) for coordinate in line.split(",")[1:]] for line in open(reference_nodes_file_basename+".inp", "r").readlines()])

    if (verbose): print "Read initial coordinates"
    X_cur = numpy.array([[float(coordinate) for coordinate in line.split(",")[1:]] for line in open(initial_nodes_file_basename+".inp", "r").readlines()])

    if (verbose): print "Start loop"
    num_iteration = 0
    while (True):
        num_iteration += 1
        if (verbose): print "num_iteration = " + str(num_iteration)

        if (verbose): print "Write current reference coordinates"
        open("tmp-nodes.inp", "w").write("\n".join([", ".join([str(num_node+1)]+[str(X) for X in X_cur[num_node]]) for num_node in range(len(X_cur))]))

        if (verbose): print 'Write input file'
        shutil.copy(input_file_basename+".inp", 'tmp-job.inp')
        for sedline in sedlines:
            os.system("sed -i "+sedline+" tmp-job.inp")

        if (verbose): print "Submit job"
        os.system("abaqus analysis background" + ((" cpus="+str(cpus)) * (cpus!=None)) + ((" gpus="+str(gpus)) * (gpus!=None)) + ((" user="+user) * (user!=None)) + " job=tmp-job")

        if (verbose): print "Wait for job to complete"
        success = wait_for_job_to_complete("tmp-job")
        if not (success):
            if (verbose): print "Job failed. Aborting"
            return False

        if (verbose): print "Extract node positions"
        os.system("abaqus python ../ABAQUS_py/extract_nodes_position.py tmp-job")

        if (verbose): print "Read node positions"
        X_def = numpy.array([[float(x) for x in line.split(",")] for line in open("tmp-job.nodes_position.dat", "r").readlines()])
        os.remove("tmp-job.nodes_position.dat")

        if (verbose): print "Save ODB"
        shutil.copy('tmp-job.odb', 'iter'+str(num_iteration)+'.odb')

        if (verbose): print "Compute residual"
        if (num_iteration > 1):
            R_old = numpy.copy(R)
        R = X_def - X_ref

        if (verbose): print "Compute error"
        error = numpy.max(numpy.abs(R))
        if (verbose): print "error = " + str(error)

        if (verbose): print "Compute relaxation"
        if (num_iteration == 1):
            alpha = 1.
        else:
            alpha = - alpha * numpy.trace(numpy.dot(R_old, numpy.transpose(R-R_old))) / numpy.trace(numpy.dot(R-R_old, numpy.transpose(R-R_old)))
        if (verbose): print "alpha = " + str(alpha)

        if (verbose): print "Exit test"
        if (error < tol):
            if (verbose): print "Converged"
            return True
        if (num_iteration == nb_max_iterations):
            if (verbose): print "Too many iterations. Aborting"
            return False

        if (verbose): print "Update current reference coordinates"
        X_cur -= alpha * R
