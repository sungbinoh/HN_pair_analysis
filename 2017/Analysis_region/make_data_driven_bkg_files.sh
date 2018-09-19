#!/bin/bash

root -l -b<<-EOF
.L make_emu_file.C+
plot()
.q
EOF

root -l -b<<-EOF
.L make_CF_file.C+
plot()
.q
EOF
