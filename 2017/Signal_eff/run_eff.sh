root -l -b<<EOF
.L signal_eff.C+
plot()
.q
EOF

scp ./pdfs/*.pdf suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/Signal_Eff
scp -r cutflow suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/Signal_Eff