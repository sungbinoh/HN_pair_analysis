:<<END
root -l -b -q "src/QuickLimit.C(0)"
root -l -b -q "src/QuickLimit.C(1)"
root -l -b -q "src/QuickLimit.C(2)"
root -l -b -q "src/QuickLimit.C(3)"
root -l -b -q "src/QuickLimit.C(4)"
root -l -b -q "src/QuickLimit.C(5)"
root -l -b -q "src/QuickLimit.C(6)"
root -l -b -q "src/QuickLimit.C(7)"
root -l -b -q "src/QuickLimit.C(8)"
root -l -b -q "src/QuickLimit.C(9)"
root -l -b -q "src/QuickLimit.C(10)"
root -l -b -q "src/QuickLimit.C(11)"
root -l -b -q "src/QuickLimit.C(12)"
root -l -b -q "src/QuickLimit.C(13)"
root -l -b -q "src/QuickLimit.C(14)"
root -l -b -q "src/QuickLimit.C(15)"
root -l -b -q "src/QuickLimit.C(16)"
root -l -b -q "src/QuickLimit.C(17)"
root -l -b -q "src/QuickLimit.C(18)"
root -l -b -q "src/QuickLimit.C(19)"
END
#root -l -b -q "src/Limit_signal_injection.C(0)"
#root -l -b -q "src/Limit_signal_injection.C(1)"
#root -l -b -q "src/Limit_signal_injection.C(3)"
#root -l -b -q "src/Limit_signal_injection.C(6)"

# 3AK8 CRs obs
#:<<END
root -l -b -q "src/QuickLimit_obs.C(0)"
root -l -b -q "src/QuickLimit_obs.C(1)"
root -l -b -q "src/QuickLimit_obs.C(2)"
root -l -b -q "src/QuickLimit_obs.C(3)"
root -l -b -q "src/QuickLimit_obs.C(4)"
root -l -b -q "src/QuickLimit_obs.C(5)"
root -l -b -q "src/QuickLimit_obs.C(6)"
root -l -b -q "src/QuickLimit_obs.C(7)"
root -l -b -q "src/QuickLimit_obs.C(8)"
root -l -b -q "src/QuickLimit_obs.C(9)"
root -l -b -q "src/QuickLimit_obs.C(10)"
root -l -b -q "src/QuickLimit_obs.C(11)"
root -l -b -q "src/QuickLimit_obs.C(12)"
root -l -b -q "src/QuickLimit_obs.C(13)"
root -l -b -q "src/QuickLimit_obs.C(14)"
root -l -b -q "src/QuickLimit_obs.C(15)"
root -l -b -q "src/QuickLimit_obs.C(16)"
root -l -b -q "src/QuickLimit_obs.C(17)"
root -l -b -q "src/QuickLimit_obs.C(18)"
root -l -b -q "src/QuickLimit_obs.C(19)"

scp -r ./plots/limit suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/
#END

# -- Comparison
root -l -b -q "src/Quick_Contour_overlap_obs.C(0)"

scp -r ./plots/limit/comparison suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/limit/


#root -l -b -q src/Comparison_Hybrid_Asym.C
#scp /data6/Users/suoh/HN_pair_analysis/plots/limit/3AK8_CRs/MuMu/Asym_vs_Hybridnew.pdf suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/limit/3AK8_CRs/MuMu
#scp /data6/Users/suoh/HN_pair_analysis/plots/limit/3AK8_CRs/EE/Asym_vs_Hybridnew.pdf suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/limit/3AK8_CRs/EE
