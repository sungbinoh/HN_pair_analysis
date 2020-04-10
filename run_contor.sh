#root -l -b -q "src/QuickLimit.C(0)"
#root -l -b -q "src/QuickLimit.C(1)"
#root -l -b -q "src/QuickLimit.C(2)"
#root -l -b -q "src/QuickLimit.C(3)"
#root -l -b -q "src/QuickLimit.C(4)"
#root -l -b -q "src/QuickLimit.C(5)"
#root -l -b -q "src/QuickLimit.C(6)"
#root -l -b -q "src/QuickLimit.C(7)"
root -l -b -q "src/QuickLimit.C(8)"
#root -l -b -q "src/QuickLimit.C(9)"
#root -l -b -q "src/QuickLimit.C(10)"
root -l -b -q "src/QuickLimit.C(11)"
#root -l -b -q "src/QuickLimit.C(12)"
#root -l -b -q "src/QuickLimit.C(13)"
#root -l -b -q "src/QuickLimit.C(14)"
#root -l -b -q "src/QuickLimit.C(15)"
#root -l -b -q "src/QuickLimit.C(16)"
#root -l -b -q "src/QuickLimit.C(17)"
#root -l -b -q "src/QuickLimit.C(18)"
#root -l -b -q "src/QuickLimit.C(19)"

#root -l -b -q "src/Limit_signal_injection.C(0)"
#root -l -b -q "src/Limit_signal_injection.C(1)"
#root -l -b -q "src/Limit_signal_injection.C(3)"
#root -l -b -q "src/Limit_signal_injection.C(6)"


#scp -r ./plots/limit suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/

root -l -b -q "src/Quick_Contour_overlap.C(0)"

scp -r ./plots/limit/comparison suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/limit/
