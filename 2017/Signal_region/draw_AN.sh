#!/bin/bash
root -l -b <<-EOF
.L plot_OS.C+
plot()
.q
EOF
echo 'next job'
root -l -b <<-EOF
.L plot_syst_all_vector.C+
plot()
.q
EOF
scp -r plots/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/CR/Top_Pt_Corr