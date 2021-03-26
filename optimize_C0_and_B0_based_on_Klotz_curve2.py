#coding=utf8

### ---------------------------------------------------- ### IMPORTS ###

import os
import time

import numpy
import sympy

import nlopt

from optimize_C0_or_Tmax_based_on_volumes import *

######################################################### EXCEPTIONS ###

class FailedComputationException(Exception): pass
class ConvergedOptimizationException(Exception): pass

############################################################### CASE ###

class case_for_C0_and_B0_optimization2:
    def __init__(self, name, V_ES, V_ED, P_dia, B1, B0_ini, B0_min, B0_max, C0_ini, C0_min, C0_max, tol=1e-3):
        self.name = name
        self.V_ES = V_ES
        self.V_ED = V_ED
        self.EF = 100*(V_ED-V_ES)/V_ED
        self.P_dia = P_dia
        self.B1 = B1
        self.B0_ini = B0_ini
        self.B0_min = B0_min
        self.B0_max = B0_max
        self.C0_ini = C0_ini
        self.C0_min = C0_min
        self.C0_max = C0_max
        self.tol = tol

###################################################### COST FUNCTION ###

def cost_function_for_C0_and_B0_optimization2(x, dx):
    assert (dx.size == 0), "dx should be empty since we're using derivative-free algorithms. Aborting."

    # num_iter_for_B0_optimization
    global num_iter_for_B0_optimization
    num_iter_for_B0_optimization += 1
    print "* num_iter_for_B0_optimization = " +  str(num_iter_for_B0_optimization) + " *"

    # B0
    B0 = x[0]
    print "B0 = " + str(B0)

    # B0 variation
    global B0_old
    if (num_iter_for_B0_optimization > 1):
        B0_rel_var = (B0-B0_old)/B0_old
        print "B0_old = " + str(B0_old)
        print "B0_rel_var = " + str(100*B0_rel_var) + " %"
    B0_old = B0

    # prepare input file
    global case_for_B0_optimization
    case = case_for_C0_optimization(case_for_B0_optimization.name,
                                    case_for_B0_optimization.V_ES,
                                    case_for_B0_optimization.V_ED,
                                    case_for_B0_optimization.P_dia,
                                    case_for_B0_optimization.B1,
                                    B0,
                                    case_for_B0_optimization.C0_ini,
                                    case_for_B0_optimization.C0_min,
                                    case_for_B0_optimization.C0_max,
                                    case_for_B0_optimization.tol)
    C0 = optimize_C0_or_Tmax_based_on_volumes([case])
    os.system("mv Error.dat Error_for_C0_optimization-Iter" + str(num_iter_for_B0_optimization) + ".dat ")
    print "C0 = " + str(C0) + " kPa"

    # extract data
    os.system("abq6132 python ../extract_fluid_cavity_volumes_and_pressures.py Heart")
    volumes = [float(volume)/1e3 for volume in open("Heart.volumes.dat", "r").read()[1:-1].split(", ")]
    pressures = [float(pressure)/0.133322 for pressure in open("Heart.pressures.dat", "r").read()[1:-1].split(", ")]
    assert (len(volumes) == len(pressures)), "Volumes and Pressures should have same length. Aborting."
    normalized_volumes = [(volume-volumes[0])/(volumes[-1]-volumes[0]) for volume in volumes]
    #print "volumes = " + str(volumes)
    #print "pressures = " + str(pressures)
    #print "normalized_volumes = " + str(normalized_volumes)

    # compute error
    global err, V1, V2, P1, P2, err0
    error = sum([err.subs(V1, V1_).subs(V2, V2_).subs(P1, P1_).subs(P2, P2_).doit() for (V1_, V2_, P1_, P2_) in zip(normalized_volumes[:-1], normalized_volumes[1:], pressures[:-1], pressures[1:])])**0.5 / err0
    error = float(error)
    print "error = " + str(error)

    # save error
    if (num_iter_for_B0_optimization == 1):
        file = open("Error_for_B0_optimization.dat", "w")
        file.write("#num_iter C0 (kPa) B0 () B1 () error ()\n")
        file.close()
    file = open("Error_for_B0_optimization.dat", "a")
    file.write(str(num_iter_for_B0_optimization) + " " + str(C0) + " " + str(B0) + " " + str(case_for_B0_optimization.B1) + " " + str(error) + "\n")
    file.close()

    # stopping criteria
    if (num_iter_for_B0_optimization > 1) and (abs(B0_rel_var) < case_for_B0_optimization.tol):
    #if (num_iter_for_B0_optimization > 1) and (error < case_for_B0_optimization.tol) and (abs(B0_rel_var) < case_for_B0_optimization.tol):
        os.system("mv Heart.odb " + case.name + "-Dia.odb")
        os.system("mv Heart.dat " + case.name + "-Dia.dat")
        os.system("rm -f Heart*")
        raise ConvergedOptimizationException

    ## stopping criteria
    #if (error < 1e-2): raise ConvergedOptimizationException

    # return
    return error

################################################################ RUN ###

def optimize_C0_and_B0_based_on_Klotz_curve2(case_):
    global case_for_B0_optimization
    case_for_B0_optimization = case_

    opt = nlopt.opt(nlopt.LN_BOBYQA, 1)

    opt.set_min_objective(cost_function_for_C0_and_B0_optimization2)
    opt.set_lower_bounds([case_for_B0_optimization.B0_min])
    opt.set_upper_bounds([case_for_B0_optimization.B0_max])

    global err, V1, V2, P1, P2, err0
    V, V1, V2, P1, P2 = sympy.symbols("V V1 V2 P1 P2")
    Pnum = P1 + (V-V1)/(V2-V1)*(P2-P1)
    Pref = case_for_B0_optimization.P_dia * V**2.76
    err  = sympy.Integral((Pnum-Pref)**2, (V,V1,V2))
    err0 = sympy.integrate(Pref**2, (V,0,1))**0.5

    os.system("rm -f Error*")
    os.system("rm -f Heart*")

    try:
        global num_iter_for_B0_optimization
        num_iter_for_B0_optimization = 0
        global B0_old
        B0_old = case_for_B0_optimization.B0_ini
        opt.optimize([case_for_B0_optimization.B0_ini])
    except ConvergedOptimizationException:
        print "ConvergedOptimizationException"
    except FailedComputationException:
        print "FailedComputationException"
