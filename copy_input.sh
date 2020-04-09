cp_2018_bkg="cp "${path_2018}'HN_pair_all_SkimTree_LRSMHighPt*.root '$FILE_PATH'/2018/Background/'
rm_signal_bkg_2018="rm "$FILE_PATH'/2018/Background/HN_pair_all_HNPairToJJJJ_*.root'
cp_2018_DATA="cp "${path_2018}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2018/DATA/'
cp_2018_signal="cp "${path_2018}'HN_pair_all_HNPairToJJJJ_*.root '$FILE_PATH'/2018/Signal/'
cp_2018_SR_ZpNN="cp "${path_SR_ZpNN}'2018/SR_ZpNN_HNPairToJJJJ_*.root '$FILE_PATH'/2018/Signal/'
echo "copying 2018 Background"
#eval $cp_2018_bkg
#eval $rm_signal_bkg_2018
echo "copying 2018 signal"
#eval $cp_2018_signal
#eval $cp_2018_SR_ZpNN
echo "copying 2018 DATA"
#eval $cp_2018_DATA

cp_2017_bkg="cp "${path_2017}'HN_pair_all_SkimTree_LRSMHighPt*.root '$FILE_PATH'/2017/Background/'
rm_signal_bkg_2017="rm "$FILE_PATH'/2017/Background/HN_pair_all_HNPairToJJJJ_*.root'
cp_2017_DATA="cp "${path_2017}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2017/DATA/'
cp_2017_signal="cp "${path_2017}'HN_pair_all_HNPairToJJJJ_*.root '$FILE_PATH'/2017/Signal/'
cp_2017_SR_ZpNN="cp "${path_SR_ZpNN}'2017/SR_ZpNN_HNPairToJJJJ_*.root '$FILE_PATH'/2017/Signal/'
echo "copying 2017 Background"
#eval $cp_2017_bkg
#eval $rm_signal_bkg_2017
echo "copying 2017 signal"
#eval $cp_2017_signal
eval $cp_2017_SR_ZpNN
echo "copying 2017 DATA"
#eval $cp_2017_DATA

cp_2016_bkg="cp "${path_2016}'HN_pair_all_SkimTree_LRSMHighPt*.root '$FILE_PATH'/2016/Background/'
rm_signal_bkg_2016="rm "$FILE_PATH'/2016/Background/HN_pair_all_HNPairToJJJJ_*.root'
cp_2016_DATA="cp "${path_2016}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2016/DATA/'
cp_2016_signal="cp "${path_2016}'HN_pair_all_HNPairToJJJJ_*.root '$FILE_PATH'/2016/Signal/'
cp_2016_SR_ZpNN="cp "${path_SR_ZpNN}'2016/SR_ZpNN_HNPairToJJJJ_*.root '$FILE_PATH'/2016/Signal/'
echo "copying 2016 Background"
#eval $cp_2016_bkg
#eval $rm_signal_bkg_2016
echo "copying 2016 signal"
#eval $cp_2016_signal
eval $cp_2016_SR_ZpNN
echo "copying 2016 DATA"
#eval $cp_2016_DATA
