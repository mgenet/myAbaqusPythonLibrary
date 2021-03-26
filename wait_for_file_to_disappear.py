#coding=utf8

### ---------------------------------------------------------------- ###

import os
import sys
import time

### ---------------------------------------------------------------- ###

def wait_for_file_to_disappear(filename, sleep=1., verbose=True):
    if (verbose): sys.stdout.write("Wait for file " + filename + " to disappear")
    while os.path.exists(filename):
        if (verbose): sys.stdout.write(".")
        time.sleep(sleep)
    if (verbose): print(".")
    time.sleep(sleep)
    return

