#root -l -b -q "src/QuickLimit.C(0)"
root -l -b -q "src/QuickLimit.C(1)"
root -l -b -q "src/QuickLimit.C(2)"
root -l -b -q "src/QuickLimit.C(3)"
#root -l -b -q "src/QuickLimit.C(4)"
root -l -b -q "src/QuickLimit.C(5)"
root -l -b -q "src/QuickLimit.C(6)"
root -l -b -q "src/QuickLimit.C(7)"

scp -r ./plots/limit suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/
