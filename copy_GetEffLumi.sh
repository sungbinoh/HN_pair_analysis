cp_2016="cp "${path_GetEffLumi}'/2016/*.root '$FILE_PATH'/2016/GetEffLumi/'
cp_2017="cp "${path_GetEffLumi}'/2017/*.root '$FILE_PATH'/2017/GetEffLumi/'
cp_2018="cp "${path_GetEffLumi}'/2018/*.root '$FILE_PATH'/2018/GetEffLumi/'

echo "copying 2016 GetEffLumi"
eval $cp_2016
echo "copying 2017 GetEffLumi"
eval $cp_2017
echo "copying 2018 GetEffLumi"
eval $cp_2018
