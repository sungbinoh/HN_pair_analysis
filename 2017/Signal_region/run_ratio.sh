#!/bin/bash
#rm HN_pair_all_VV.root
#hadd HN_pair_all_VV.root HN_pair_all_ZZTo4L_powheg.root HN_pair_all_ZZTo2L2Q.root HN_pair_all_WZTo3LNu.root HN_pair_all_WZTo2L2Q.root HN_pair_all_WWTo2L2Nu_powheg.root HN_pair_all_ZZZ.root HN_pair_all_WZZ.root HN_pair_all_WWZ.root HN_pair_all_WWW.root HN_pair_all_SingleTop_tch_top.root HN_pair_all_SingleTop_tch_antitop.root HN_pair_all_SingleTop_tW_top.root HN_pair_all_SingleTop_tW_antitop.root HN_pair_all_SingleTop_sch.root HN_pair_all_TTG.root HN_pair_all_WJets_MG.root HN_pair_all_ttZ.root HN_pair_all_ttW.root

root -l -b <<-EOF
.L plot_prompt_level.C+
plot()
.q
EOF
echo 'next job'
scp -r plots/ratio suoh@147.47.242.44:/home/suoh/Dropbox/HN_pair/2017/CR/Syst
