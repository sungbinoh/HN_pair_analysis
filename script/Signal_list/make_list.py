import os
import sys
import commands
import datetime
import subprocess
import time


year = sys.argv[1]

#year = "2016"

a = 0
f = open("./Zprime_to_NN_" + year + ".txt", 'r')

for line in f: # -- loop over mass points 
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    print mass_point
    channel = mass_point.split("_")[1][0:2]
    if "Mu" in channel:
        channel = channel + channel
    print channel
    mZp = mass_point.split("_")[2][6:]
    print mZp
    mN = mass_point.split("_")[3][1:]
    print mN
    list_name = "HNPairToJJJJ_" + channel + "_ZP" + mZp + "_N" + mN + "_WR5000.txt"
    print list_name
    f_current_list = open("./list_" + year + "/" + list_name, 'w')
    
    current_target_dir = "/gv0/DATA/SKFlat/Run2Legacy_v4/" + year + "/PrivateMC/ZpNN/" + mass_point
    print current_target_dir
    
    file_list = commands.getstatusoutput("ls " + current_target_dir)
    #print file_list
    file_list_split = file_list[1].split("\n")
    size_file_list_split = len(file_list_split)
    #print file_list_split
    #print size_file_list_split
    for i in range(0, size_file_list_split): 
        f_current_list.write(current_target_dir + "/" + file_list_split[i] + "\n")
        


    f_current_list.close()
    #abort for test
    #a = a + 1
    #if a > 5:
    #    break
f.close()
