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

cp -r emu_ratio/* /home/suoh/Dropbox/HN_pair/2016/Emu_ratio
