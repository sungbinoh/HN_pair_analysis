#!/bin/bash

root -l -b <<EOF
.L plot_for_AN.C+
plot()
.q
EOF

scp -r ./plots/* suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/CR