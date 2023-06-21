#!/bin/bash
rm HN_pair_all_VV.root
hadd HN_pair_all_VV.root HN_pair_all_ZZTo4L_powheg.root HN_pair_all_ZZTo2L2Q.root HN_pair_all_WZTo3LNu.root HN_pair_all_WZTo2L2Q.root HN_pair_all_WWTo2L2Nu_powheg.root

root -l -b <<-EOF
.L plot_OS.C+
plot()
.q
EOF
echo 'next job'
#root -l -b <<-EOF
#.L plot_SS.C+
#plot()
#.q
#EOF
scp -r plots/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/CR
