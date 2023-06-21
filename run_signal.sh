#root -l -b -q "src/Prefire_SignalEff_Effect.C(2016)"
#root -l -b -q "src/Prefire_SignalEff_Effect.C(2017)"

#root -l -b -q "src/Signal_private_VS_official.C(2016)"
#root -l -b -q "src/Signal_private_VS_official.C(2017)"
#root -l -b -q "src/Signal_private_VS_official.C(2018)"

#root -l -b -q "src/QuickSignal.C(2016)"

root -l -b -q src/HN_matching.C

scp -r ./plots/Signal/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/Signal
