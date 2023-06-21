setup_root='sr ~/setup_root.sh'
rm_VV='rm HN_pair_all_VV.root'
rm_old_DY='rm HN_pair_all_DYJets.root'
SingleMuon='mv HN_pair_all_data_SingleMuon_cat_v8-0-8.root HN_pair_all_SingleMuon.root'
DoubleEG='mv HN_pair_all_data_DoubleEG_cat_v8-0-8.root HN_pair_all_DoubleEG.root'
rm_TT='rm HN_pair_all_SKTT_powheg_dilep_cat_v8-0-8.root '
rm_old_TT='rm HN_pair_all_TTLL_powheg.root'
TTLL='mv HN_pair_all_SKTTLL_powheg_dilep_cat_v8-0-8.root HN_pair_all_TTLL_powheg.root'
#DY='hadd HN_pair_all_DYJets.root HN_pair_all_SKDY10to50_dilep_cat_v8-0-8.root HN_pair_all_SKDY50plus_dilep_cat_v8-0-8.root'
DY='hadd HN_pair_all_DYJets.root HN_pair_all_SKDY_pt_*'
rm_DY='rm HN_pair_all_SKDY*cat_v8-0-8.root'
Wjet='mv HN_pair_all_SKWJets_dilep_cat_v8-0-8.root HN_pair_all_WJets_MG.root'
other='hadd HN_pair_all_VV.root HN_pair_all_SK*.root'
rm_other='rm *v8-0-8.root'

eval $rm_VV
eval $rm_old_DY

eval $setup_root

echo $SingleMuon
eval $SingleMuon

echo $DoubleEG
eval $DoubleEG

eval $rm_old_TT
eval $rm_TT
echo $TTLL
eval $TTLL

echo $DY
eval $DY
echo $rm_DY
eval $rm_DY

echo $Wjet
eval $Wjet

echo $other
eval $other
echo $rm_other
eval $rm_other