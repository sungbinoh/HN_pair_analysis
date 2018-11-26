setup_root='sr ~/setup_root.sh'
SingleMuon='mv HN_pair_all_data_SingleMuon_cat_v8-0-8.root HN_pair_all_SingleMuon.root'
DoubleEG='mv HN_pair_all_data_DoubleEG_cat_v8-0-8.root HN_pair_all_DoubleEG.root'
TTLL='mv HN_pair_all_SKTT_powheg_dilep_cat_v8-0-8.root HN_pair_all_TTLL_powheg.root'
DY='hadd HN_pair_all_DYJets.root HN_pair_all_SKDY10to50_dilep_cat_v8-0-8.root HN_pair_all_SKDY50plus_dilep_cat_v8-0-8.root'
rm_DY='rm HN_pair_all_SKDY*cat_v8-0-8.root'
Wjet='mv HN_pair_all_SKWJets_dilep_cat_v8-0-8.root HN_pair_all_WJets_MG.root'
other='hadd HN_pair_all_VV.root HN_pair_all_SK*.root'
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

echo $Wjet
eval $Wjet

echo $other
eval $other
echo $rm_other
eval $rm_other