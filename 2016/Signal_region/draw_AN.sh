#!/bin/bash
root -l -b <<-EOF
.L plot_OS.C+
plot()
.q
EOF
#root -l -b <<-EOF
#.L plot_mZp.C+
#plot()
#.q
#EOF
#root -l -b <<-EOF
#.L plot_syst_all.C+
#plot()
#.q
#EOF
#root -l -b <<-EOF
#.L plot_syst_all_mZp.C+
#plot()
#.q
#EOF
scp -r plots/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2016/CR/Syst
