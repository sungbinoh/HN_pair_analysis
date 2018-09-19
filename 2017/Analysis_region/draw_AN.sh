#!/bin/bash

root -l -b <<-EOF
.L plot_OS.C+
plot()
.q
EOF
echo 'next job'
#root -l -b <<-EOF
#.L plot_SS_for_AN_emu_method_DY_amcnlo.C+
#plot()
#.q
#EOF
scp -r plots/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/CR
