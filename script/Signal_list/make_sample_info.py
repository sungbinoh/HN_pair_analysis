import os
import sys
import commands
import datetime
import subprocess
import time


year = sys.argv[1]

#year = "2016"

a = 0
#f = open("./Zprime_to_NN_" + year + ".txt", 'r')
f = open("./Zprime_to_NN_2016.txt", 'r')

for line in f: # -- loop over mass points 
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    print mass_point
    channel = mass_point.split("_")[1][0:2]
    print channel
    if "Mu" in channel:
         channel = channel + channel
    mZp = mass_point.split("_")[2][6:]
    print mZp
    mN = mass_point.split("_")[3][1:]
    print mN
    list_name = "HNPairToJJJJ_" + channel + "_ZP" + mZp + "_N" + mN + "_WR5000.txt"
    print list_name
    f_current_list = open("./CommonSampleInfo/" + year + "/" + list_name, 'w')
    
    f_current_list.write("# alias PD xsec nmc sumw\n")
    f_current_list.write("HNPairToJJJJ_" + channel + "_ZP" + mZp + "_N" + mN + "_WR5000" + "\t" + mass_point + "\t0.1\t50000\t50000")
    
    f_current_list.close()
    #abort for test
    #a = a + 1
    #if a > 5:
    #    break
f.close()
