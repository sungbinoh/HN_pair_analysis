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
cp -r plots_emu_method/* /home/suoh/Dropbox/HN_pair/2016/CR
