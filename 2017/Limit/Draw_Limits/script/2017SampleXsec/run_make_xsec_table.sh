#!/bin/bash

root -l -b<<EOF
.L make_xsec_table.C+
plot()
.q
EOF

