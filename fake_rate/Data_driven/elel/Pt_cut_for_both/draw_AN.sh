#!/bin/bash

root -l -b <<-EOF
.L plot_for_AN.C+
plot()
.q
EOF
echo 'next job'
#root -l -b <<-EOF
#.L plot_for_AN_log.C+
#plot()
#.q
#EOF
cp pdfs_HLT_CR/*.pdf /home/suoh/Dropbox/HN_pair/2016/fake_rate/Electron/Pt_cut_for_both/kinematics
