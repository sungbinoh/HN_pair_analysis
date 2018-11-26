#!/bin/bash

root -l -b<<EOF
.L make_nobin_OS.C+
plot()
.q
EOF

#root -l -b<<EOF
#.L make_nobin_SS.C+
#plot()
#.q
#EOF
