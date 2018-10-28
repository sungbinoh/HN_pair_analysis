#!/bin/bash

root -l -b <<-EOF
.L plot_for_AN.C+
plot()
.q
EOF
echo 'next job'
root -l -b <<-EOF
.L plot_for_AN_log.C+
plot()
.q
EOF
scp pdfs_HLT_CR/*.pdf suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2016/fake_rate/Electron/kinematics
