#!/bin/bash

root -l -b <<-EOF
.L Get_FR.C+
plot()
.q
EOF
cp ./Fake_Rate/*.pdf /home/suoh/Dropbox/HN_pair/2016/fake_rate/2D_rate

