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
scp ./Fake_Rate/*.pdf suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2016/fake_rate/Electron/2D_rate
