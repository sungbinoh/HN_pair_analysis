#!/bin/bash

root -l -b <<-EOF
.L Get_FR.C+
plot()
.q
EOF
root -l -b <<-EOF
.L Make_FR_plot_file.C+
plot()
.q
EOF
root -l -b <<-EOF
.L Draw_FR_1D.C+
plot()
.q
EOF
cp ./Fake_Rate/*.pdf /home/suoh/Dropbox/HN_pair/2016/fake_rate/Muon/v8-0-8/2D_rate
