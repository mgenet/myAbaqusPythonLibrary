#coding=utf8

### ---------------------------------------------------- ### IMPORTS ###

import os
import time

import numpy
import sympy

import nlopt

from .wait_for_job_to_complete import *

######################################################### EXCEPTIONS ###

class FailedComputationException(Exception): pass
class ConvergedOptimizationException(Exception): pass

###################################################### COST FUNCTION ###

def cost_function_for_C0_and_B0_optimization(params, dx):
    assert (dx.size == 0), "dx should be empty since we're using derivative-free algorithms. Aborting."

    # num_iter
    global num_iter
    num_iter += 1
    print("* num_iter = " +  str(num_iter) + " *")

    # params
    C0 = params[0]
    B0 = params[1]
    if (len(params) == 3):
        D0 = params[2]
    else:
        D0 = 0.40
    print("C0 = " + str(C0) + " kPa")
    print("B0 = " + str(B0))
    if (len(params) == 3): print("D0 = " + str(D0))

    # params variation
    global params_old
    if (num_iter > 1):
        print("C0_old = " + str(params_old[0]) + " kPa")
        print("B0_old = " + str(params_old[1]))
        if (len(params) == 3): print("D0_old = " + str(params_old[2]))
        params_rel_var = [(param-param_old)/param_old for (param, param_old) in zip(params, params_old)]
        print("C0_rel_var = " + str(100*params_rel_var[0]) + " %")
        print("B0_rel_var = " + str(100*params_rel_var[1]) + " %")
        if (len(params) == 3): print("D0_rel_var = " + str(100*params_rel_var[2]) + " %")
    params_old = [param_old for param_old in params]

    # prepare input file
    global name
    os.system("rm -f Heart*")
    os.system("cp ../Heart_input.inp Heart.inp")
    os.system("sed -i 's/<<Name>>/" + name + "/g' Heart.inp")
    os.system("sed -i 's/<<C0>>/" + str(C0) + "/g' Heart.inp")
    os.system("sed -i 's/<<B01111>>/" + str(B0) + "/g' Heart.inp")
    os.system("sed -i 's/<<B02222>>/" + str(B0*D0) + "/g' Heart.inp")
    os.system("sed -i 's/<<B03333>>/" + str(B0*D0) + "/g' Heart.inp")
    os.system("sed -i 's/<<B01212>>/" + str(B0*(1.+D0)/4) + "/g' Heart.inp")
    os.system("sed -i 's/<<B01313>>/" + str(B0*(1.+D0)/4) + "/g' Heart.inp")
    os.system("sed -i 's/<<B02323>>/" + str(B0*D0/2) + "/g' Heart.inp")

    # run job
    os.system("abq6132 user=../usub.o job=Heart")

    # wait for job to complete
    success = wait_for_job_to_complete("Heart")

    # compute error
    if (success):
        # extract data
        os.system("abq6132 python ../extract_fluid_cavity_volumes_and_pressures.py Heart")
        volumes = [float(volume)/1e3 for volume in open("Heart.volumes.dat", "r").read()[1:-1].split(", ")]
        ejection_fraction = 100 * (volumes[-1] - volumes[0]) / volumes[-1]
        print("ejection fraction = " + str(ejection_fraction))
        pressures = [float(pressure)/0.133322 for pressure in open("Heart.pressures.dat", "r").read()[1:-1].split(", ")]
        assert (len(volumes) == len(pressures)), "Volumes and Pressures should have same length. Aborting."
        normalized_volumes = [(volume-volumes[0])/(volumes[-1]-volumes[0]) for volume in volumes]
        #print "volumes = " + str(volumes)
        #print "pressures = " + str(pressures)
        #print "normalized_volumes = " + str(normalized_volumes)

        # save data
        file = open("Iter" + str(num_iter) + ".dat", "w")
        file.write("# volume normalized_volume pressure\n")
        for volume, normalized_volume, pressure in zip(volumes, normalized_volumes, pressures):
            file.write(str(volume) + " " + str(normalized_volume) + " " + str(pressure) + "\n")
        file.close()

        # print data
        file = open("Iter.gnu", "w")
        file.write('''\
set term pdf enhanced
set output "Iter.pdf"
set xlabel "normalized volume ()"
set ylabel "pressure (mmHg)"
''')
        for k_iter in range(1,num_iter+1):
            file.write('''\
''' + ((k_iter==1) * '''plot ''') + ((k_iter>1) * '''     ''') + '''"Iter''' + str(k_iter) + '''.dat" using ($2):($3) with lines linewidth 3 title "Iter''' + str(k_iter) + '''",\\
''')
        file.write('''\
     10*x**2.76 with lines linecolor "black" linewidth 5 title "[Klotz et al., 2006, Am. J. Physiol. Heart. Circ. Physiol.]"
''')
        file.close()
        os.system("gnuplot Iter.gnu")

        # compute error
        global err, V1, V2, P1, P2, err0
        error_PV = sum([err.subs(V1, V1_).subs(V2, V2_).subs(P1, P1_).subs(P2, P2_).doit() for (V1_, V2_, P1_, P2_) in zip(normalized_volumes[:-1], normalized_volumes[1:], pressures[:-1], pressures[1:])])**0.5 / err0
        error_PV = float(error_PV)
        print("error_PV = " + str(error_PV))

        global target
        error_EF = abs(ejection_fraction-target)/target
        print("error_EF = " + str(error_EF))

        error = error_PV + error_EF
        print("error = " + str(error))
    else:
        error_PV = 1.
        ejection_fraction = 0.
        error_EF = 1.
        error = 2.

    # save error
    if (num_iter == 1):
        file = open("Error.dat", "w")
        file.write("#num_iter C0 (kPa) B0 ()" + (len(params) == 3)*(" D0 ()") + " ejection fraction (%) error_PV () error_EF () error ()\n")
        file.close()
    file = open("Error.dat", "a")
    file.write(str(num_iter) + " " + str(C0) + " " + str(B0) + (len(params) == 3)*(" " + str(D0)) + " " + str(ejection_fraction) + " " + str(error_PV) + " " + str(error_EF) + " " + str(error) + "\n")
    file.close()

    # stopping criteria
    if (num_iter > 1) and (error < 1e-3) and (max(abs(numpy.array(params_rel_var))) < 1e-3):
        os.system("rm -f Heart*")
        raise ConvergedOptimizationException

    ## stopping criteria
    #if (error < 1e-2): raise ConvergedOptimizationException

    # return
    return error

################################################################ RUN ###

def optimize_C0_and_B0_based_on_Klotz_curve(name_, target_, params, algo=nlopt.GN_DIRECT):
    opt = nlopt.opt(algo, len(params))

    global name
    name = name_

    global target
    target = target_

    global cost_function_for_C0_and_B0_optimization
    opt.set_min_objective(cost_function_for_C0_and_B0_optimization)
    opt.set_lower_bounds([param[1] for param in params])
    opt.set_upper_bounds([param[2] for param in params])

    global err, V1, V2, P1, P2, err0
    V, V1, V2, P1, P2 = sympy.symbols("V V1 V2 P1 P2")
    Pnum = P1 + (V-V1)/(V2-V1)*(P2-P1)
    Pref = 10.*V**2.76
    err  = sympy.Integral((Pnum-Pref)**2, (V,V1,V2))
    err0 = sympy.integrate(Pref**2, (V,0,1))**0.5

    os.system("rm -f Error*")
    os.system("rm -f Heart*")
    os.system("rm -f Iter*")

    try:
        global num_iter
        global params_old
        num_iter = 0
        params_old = [0. for param in params]
        opt.optimize([param[0] for param in params])
    except ConvergedOptimizationException:
        print("ConvergedOptimizationException")
    except FailedComputationException:
        print("FailedComputationException")
