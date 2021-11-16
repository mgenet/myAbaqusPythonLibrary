#coding=utf8

### ---------------------------------------------------------------- ###

from .wait_for_file_to_appear           import *
from .wait_for_file_to_disappear        import *
from .wait_for_string_to_appear_in_file import *

### ---------------------------------------------------------------- ###

def wait_for_job_to_complete(job_name, sleep=1., verbose=True):

    wait_for_file_to_appear(job_name + ".log", sleep, verbose)
    wait_for_string_to_appear_in_file((("End Abaqus/Standard Analysis", job_name + ".log"),\
                                       ("Abaqus/Analysis exited with errors", job_name + ".log"),\
                                       ("slurmstepd: error:", job_name + ".log")), sleep, verbose)
    if ("End Abaqus/Standard Analysis" in open(job_name + ".log").read()):
        success = True
    else:
        success = False
    wait_for_file_to_disappear(job_name + ".lck", sleep, verbose)
    return success
