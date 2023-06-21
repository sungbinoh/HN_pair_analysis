#!/bin/bash

root -l -b -q "src/QuickLimit.C(0)" ## MuMu, shape
#root -l -b -q "src/QuickLimit.C(1)" ## MuMU, count
root -l -b -q "src/QuickLimit.C(2)" ## ElEl, shape
#root -l -b -q "src/QuickLimit.C(3)" ## ElEl, count

scp -r plots/v946p1_3/QuickLimit/ suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017