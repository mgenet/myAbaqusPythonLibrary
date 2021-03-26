#coding=utf8

### ---------------------------------------------------- ### IMPORTS ###

import os
import time

import numpy
import nlopt

from wait_for_job_to_complete import *

######################################################### EXCEPTIONS ###

class FailedComputationException(Exception): pass
class ConvergedOptimizationException(Exception): pass

############################################################### CASE ###

class case_for_C0_optimization:
    def __init__(self, name, V_ES, V_ED, P_dia, B0, B1, B2, C0_ini, C0_min, C0_max, tol=1e-3):
        assert ("Cut" in name) or ("Full" in name), "name must contain Cut or Full. Aborting."
        self.name = name
        self.dia_or_sys = "dia"
        self.V_ES = V_ES
        self.V_ED = V_ED
        self.EF = 100*(V_ED-V_ES)/V_ED
        self.target = P_dia
        self.B0 = B0
        self.B1 = B1
        self.B2 = B2
        self.param_ini = C0_ini
        self.param_min = C0_min
        self.param_max = C0_max
        self.tol = tol

    def disp(self):
        print "*** name = " + self.name
        print "*** dia_or_sys = dia ***"
        print "*** V_ES = " + str(self.V_ES) + " mm³ ***"
        print "*** V_ED = " + str(self.V_ED) + " mm³ ***"
        print "*** EF = " + str(self.EF) + " % ***"
        print "*** P_dia = " + str(self.target) + " mmHg ***"
        print "*** B0 = " + str(self.B0) + " ***"
        print "*** B1 = " + str(self.B1) + " ***"
        print "*** B2 = " + str(self.B2) + " ***"
        print "*** C0_ini = " + str(self.param_ini) + " kPa ***"
        print "*** C0_min = " + str(self.param_min) + " kPa ***"
        print "*** C0_max = " + str(self.param_max) + " kPa ***"

    def param_name(self):
        return "C0"

class case_for_Tmax_optimization:
    def __init__(self, name, EF, C0, B0, B1, B2, P_sys, Tmax_ini, Tmax_min, Tmax_max, tol=1e-3):
        assert ("Cut" in name) or ("Full" in name), "name must contain Cut or Full. Aborting."
        self.name = name
        self.dia_or_sys = "sys"
        self.EF = EF
        self.C0 = C0
        self.B0 = B0
        self.B1 = B1
        self.B2 = B2
        self.target = P_sys
        self.param_ini = Tmax_ini
        self.param_min = Tmax_min
        self.param_max = Tmax_max
        self.tol = tol

    def disp(self):
        print "*** name = " + self.name
        print "*** dia_or_sys = sys ***"
        print "*** EF = " + str(self.EF) + " % ***"
        print "*** C0 = " + str(self.C0) + " kPa ***"
        print "*** B0 = " + str(self.B0) + " ***"
        print "*** B1 = " + str(self.B1) + " ***"
        print "*** B2 = " + str(self.B2) + " ***"
        print "*** P_sys = " + str(self.target) + " mmHg ***"
        print "*** Tmax_ini = " + str(self.param_ini) + " kPa ***"
        print "*** Tmax_min = " + str(self.param_min) + " kPa ***"
        print "*** Tmax_max = " + str(self.param_max) + " kPa ***"

    def param_name(self):
        return "Tmax"

###################################################### COST FUNCTION ###

def cost_function_error_for_C0_or_Tmax_optimization(x, dx):
    assert (dx.size == 0), "dx should be empty since we're using derivative-free algorithms. Aborting."

    # sol
    global sol
    sol = x[0]

    # num_iter
    global num_iter
    num_iter += 1
    print "* num_iter = " +  str(num_iter) + " *"

    # param
    param = x[0]
    global case
    print case.param_name() + " = " + str(param) + " kPa"

    # parameter variation
    global param_old
    if (num_iter > 1):
        print case.param_name() + "_old = " + str(param_old) + " kPa"
        if (param_old > 0):
            param_rel_var = (param-param_old)/param_old
        else:
            param_rel_var = 1.
        print case.param_name() + "_rel_var = " + str(100*param_rel_var) + " %"

    # run job
    os.system("rm -f Heart.*")
    os.system("cp ../Heart_input.inp Heart.inp")
    os.system("sed -i 's/<<Name>>/" + case.name + "/g' Heart.inp")
    if (case.dia_or_sys == "dia"):
        os.system("sed -i 's/<<C0>>/" + str(param) + "/' Heart.inp")
        os.system("sed -i 's/<<B01111>>/" + str(case.B0) + "/g' Heart.inp")
        os.system("sed -i 's/<<B02222>>/" + str(case.B1) + "/g' Heart.inp")
        os.system("sed -i 's/<<B03333>>/" + str(case.B2) + "/g' Heart.inp")
        os.system("sed -i 's/<<B01212>>/" + str((case.B0+case.B1)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<B01313>>/" + str((case.B0+case.B2)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<B02323>>/" + str((case.B1+case.B2)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<Tmax>>/" + str(0) + "/' Heart.inp")
        os.system("sed -i 's/<<Dt>>/" + str(350e-5) + "/' Heart.inp")
        os.system("sed -i 's/<<Mdot>>/" + str(1.060e-6*float(case.V_ED-case.V_ES)/350e-3) + "/' Heart.inp")
    elif (case.dia_or_sys == "sys"):
        os.system("sed -i 's/<<C0>>/" + str(case.C0) + "/' Heart.inp")
        os.system("sed -i 's/<<B01111>>/" + str(case.B0) + "/g' Heart.inp")
        os.system("sed -i 's/<<B02222>>/" + str(case.B1) + "/g' Heart.inp")
        os.system("sed -i 's/<<B03333>>/" + str(case.B2) + "/g' Heart.inp")
        os.system("sed -i 's/<<B01212>>/" + str((case.B0+case.B1)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<B01313>>/" + str((case.B0+case.B2)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<B02323>>/" + str((case.B1+case.B2)/4) + "/g' Heart.inp")
        os.system("sed -i 's/<<Tmax>>/" + str(param) + "/' Heart.inp")
        os.system("sed -i 's/<<Dt>>/" + str(350e-6) + "/' Heart.inp")
        os.system("sed -i 's/<<Mdot>>/" + str(0.) + "/' Heart.inp")
    if ("Cut" in case.name):
        os.system("sed -i 's/** N-BASEND/N-BASEND/' Heart.inp")
        os.system("sed -i 's/** N-BASMYO/N-BASMYO/' Heart.inp")
        os.system("sed -i 's/** N-BASEPI/N-BASEPI/' Heart.inp")
    elif ("Full" in case.name):
        os.system("sed -i 's/** N-BASEDG/N-BASEDG/' Heart.inp")
    os.system("abq6132 user=../usub.o job=Heart")

    # wait for job to complete
    success = wait_for_job_to_complete("Heart")
    if not success: raise FailedComputationException

    # extract data
    os.system("abq6132 python ../extract_fluid_cavity_final_pressure.py Heart")
    P = float(open('Heart.pressure.dat', 'r').read())/0.133322
    print "P = " + str(P) + " mmHg"

    # compute error
    error = (P-case.target)/case.target
    print "P_err = " + str(100*error) + " %"

    # save error
    if (num_iter == 1):
        file = open("Error.dat", "w")
        file.write("#num_iter param (kPa) error ()\n")
        file.close()
    file = open("Error.dat", "a")
    file.write(str(num_iter) + " " + str(param) + " " + str(error) + "\n")
    file.close()

    # stopping criteria
    if (num_iter > 1) and (abs(error) < case.tol) and (abs(param_rel_var) < case.tol):
        #if (case.dia_or_sys == "dia"):
            #os.system("mv Heart.odb " + case.name + "-Dia-EF=" + str(int(round(case.EF))) + "%-Pdia=" + str(int(round(case.target))) + "mmHg.odb")
            #os.system("mv Heart.dat " + case.name + "-Dia-EF=" + str(int(round(case.EF))) + "%-Pdia=" + str(int(round(case.target))) + "mmHg.dat")
        #elif (case.dia_or_sys == "sys"):
            #os.system("mv Heart.odb " + case.name + "-Sys-Psys=" + str(int(round(case.target))) + "mmHg.odb")
            #os.system("mv Heart.dat " + case.name + "-Sys-Psys=" + str(int(round(case.target))) + "mmHg.dat")
        #os.system("rm -f Heart.*")
        raise ConvergedOptimizationException
    param_old = param

    ## stopping criteria
    #if (error < 1e-2): raise ConvergedOptimizationException

    # return
    return abs(error)

################################################################ RUN ###

def optimize_C0_or_Tmax_based_on_volumes(cases):
    global case
    for case in cases:
        case.disp()

        #opt = nlopt.opt(nlopt.LN_COBYLA, 1)
        opt = nlopt.opt(nlopt.LN_BOBYQA, 1)

        opt.set_min_objective(cost_function_error_for_C0_or_Tmax_optimization)
        opt.set_lower_bounds([case.param_min])
        opt.set_upper_bounds([case.param_max])

        os.system("rm -f Heart*")

        try:
            global num_iter
            global param_old
            num_iter = 0
            param_old = 0.
            opt.optimize([case.param_ini])
        except ConvergedOptimizationException:
            print "ConvergedOptimizationException"
            global sol
            return sol
        except FailedComputationException:
            print "FailedComputationException"
            raise FailedComputationException
