#root -l -b -q "src/QuickPlot.C(2018)"
#root -l -b -q "src/QuickPlot.C(2017)"
#root -l -b -q "src/QuickPlot.C(2016)"
#root -l -b -q "src/QuickPlotFullRun2.C(2019)"
#root -l -b -q "src/QuickPlotFullRun2_signal.C(2019)"

#root -l -b -q "src/QuickPlot_signal.C(2018)"
#root -l -b -q "src/QuickPlot_signal.C(2017)"
#root -l -b -q "src/QuickPlot_signal.C(2016)"

#root -l -b -q "src/QuickPlot_CR_signal.C(2018)"
#root -l -b -q "src/QuickPlot_CR_signal.C(2017)"
#root -l -b -q "src/QuickPlot_CR_signal.C(2016)"

#root -l -b -q "src/QuickPlot_signal.C(2018)"

#root -l -b -q "src/Plot_Post_fit.C(2016)"
#root -l -b -q "src/Plot_Post_fit.C(2017)"
#root -l -b -q "src/Plot_Post_fit.C(2018)"

root -l -b -q "src/Plot_Post_fit_FullRun2.C(2019)"
root -l -b -q "src/Plot_Post_fit_FullRun2_signal.C(2019)"

#root -l -b -q "src/Plot_signal_syst_shape.C(2016)"
#root -l -b -q "src/ttbar_comparison.C(2019)"
#root -l -b -q "src/Compare_2017_Data.C(2017)"

#scp -r ./plots/2018/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2018/
#scp -r ./plots/2017/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/
#scp -r ./plots/2016/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2016/
scp -r ./plots/2019/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2019/
#scp -r ./plots/Signal/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/Signal/

#root -l -b -q src/Plot_RateParam.C
#scp ./plots/limit/Rate_param/*.pdf suoh@147.47.242.44://home/suoh/Dropbox/HN_pair/limit/Rate_param/
