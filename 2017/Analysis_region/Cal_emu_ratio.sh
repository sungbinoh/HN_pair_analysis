#!/bin/bash

root -l -b <<-EOF 
.L Get_emu_ratio.C+ 
plot()
.q
EOF
root -l -b <<-EOF
.L Get_emu_ratio_poly1.C+
plot()
.q
EOF

scp -r emu_ratio/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/EMu_ratio
