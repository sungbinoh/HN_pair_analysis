setup_root='sr ~/setup_root.sh'
SingleMuon='mv HN_pair_all_data_SingleMuon_cat_v8-0-8.root HN_pair_all_SingleMuon.root'
DoubleEG='mv HN_pair_all_data_DoubleEG_cat_v8-0-8.root HN_pair_all_DoubleEG.root'
TTLL='mv HN_pair_all_TTLL_powheg_cat_v8-0-8.root HN_pair_all_TTLL_powheg.root'
DY='hadd HN_pair_all_DYJets.root HN_pair_all_DY50plus_cat_v8-0-8.root HN_pair_all_DY10to50_cat_v8-0-8.root'
rm_DY='rm HN_pair_all_DY*cat_v8-0-8.root'
other='hadd HN_pair_all_VV.root HN_pair_all_ZZTo4L_powheg_cat_v8-0-8.root HN_pair_all_WZZ_cat_v8-0-8.root HN_pair_all_WZTo3LNu_powheg_cat_v8-0-8.root HN_pair_all_WWZ_cat_v8-0-8.root HN_pair_all_WWW_cat_v8-0-8.root HN_pair_all_WWTo2L2Nu_cat_v8-0-8.root HN_pair_all_WJets_cat_v8-0-8.root HN_pair_all_ttZ_cat_v8-0-8.root HN_pair_all_SingleTop_tW_cat_v8-0-8.root HN_pair_all_SingleTop_t_cat_v8-0-8.root HN_pair_all_SingleTop_s_cat_v8-0-8.root HN_pair_all_SingleTbar_tW_cat_v8-0-8.root HN_pair_all_SingleTbar_t_cat_v8-0-8.root'
rm_other='rm *v8-0-8.root'

eval $setup_root

echo $SingleMuon
eval $SingleMuon

echo $DoubleEG
eval $DoubleEG

echo $TTLL
eval $TTLL

echo $DY
eval $DY
echo $rm_DY
eval $rm_DY

echo $other
eval $other
echo $rm_other
eval $rm_other