#!/bin/bash

root -l -b <<-EOF
.L plot_for_AN.C+
plot()
.q
EOF
cp pdfs_HLT_CR/*.pdf /home/suoh/Dropbox/HN_pair/2016/fake_rate/kinematics
