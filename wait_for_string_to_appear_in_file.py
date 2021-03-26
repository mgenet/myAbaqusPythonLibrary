#coding=utf8

### ---------------------------------------------------------------- ###

import sys
import time

### ---------------------------------------------------------------- ###

def wait_for_string_to_appear_in_file(couples, sleep, verbose=True):
    if (verbose): sys.stdout.write("Wait for strings to appear in files " + str(couples))
    while True:
        if (verbose): sys.stdout.write(".")
        for (string,filename) in couples:
            if (string in open(filename).read()):
                if (verbose): print(".")
                time.sleep(sleep)
                return
        time.sleep(sleep)
