cp_2018_bkg="cp "${path_2018}'*.root '$FILE_PATH'/2018/Background/'
cp_2018_DATA="cp "${path_2018}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2018/DATA/'
echo "copying 2018 Background"
eval $cp_2018
echo "copying 2018 DATA"
eval $cp_2018_DATA

cp_2017_bkg="cp "${path_2017}'*.root '$FILE_PATH'/2017/Background/'
cp_2017_DATA="cp "${path_2017}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2017/DATA/'
echo "copying 2017 Background"
eval $cp_2017
echo "copying 2017 DATA"
eval $cp_2017_DATA

cp_2016_bkg="cp "${path_2016}'*.root '$FILE_PATH'/2016/Background/'
cp_2016_DATA="cp "${path_2016}'/DATA/HN_pair_all_SkimTree_LRSMHighPt_data*.root '$FILE_PATH'/2016/DATA/'
echo "copying 2016 Background"
eval $cp_2016
echo "copying 2016 DATA"
eval $cp_2016_DATA
