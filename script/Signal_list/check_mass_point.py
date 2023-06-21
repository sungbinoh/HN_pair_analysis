import os
import sys
import commands
import datetime
import subprocess
import time


#year = sys.argv[1]
#year = "2016"

f_2017 = open("./Zprime_to_NN_2017.txt", 'r')
f_2018 = open("./Zprime_to_NN_2018.txt", 'r')
all_mass_points_2017 = [""]
all_mass_points_2018 = [""]
for line in f_2017: # -- loop over mass points 
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    channel = mass_point.split("_")[1][0:2]
    if "Mu" in channel:
         channel = channel + channel
    mZp = mass_point.split("_")[2][6:]
    mN = mass_point.split("_")[3][1:]
    current_mass_point = channel + "_" + mZp + "_" + mN
    all_mass_points_2017.append(current_mass_point)
    
for line in f_2018: # -- loop over mass points 
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    channel = mass_point.split("_")[1][0:2]
    if "Mu" in channel:
         channel = channel + channel
    mZp = mass_point.split("_")[2][6:]
    mN = mass_point.split("_")[3][1:]
    current_mass_point = channel + "_" + mZp + "_" + mN
    all_mass_points_2018.append(current_mass_point)

f_2016 = open("./Zprime_to_NN_2016.txt", 'r')
print "=========== 2017 missed mass points ==========="
for line in f_2016:
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    channel = mass_point.split("_")[1][0:2]
    if "Mu" in channel:
         channel = channel + channel
    mZp = mass_point.split("_")[2][6:]
    mN = mass_point.split("_")[3][1:]
    current_mass_point = channel + "_" + mZp + "_" + mN
    if current_mass_point not in all_mass_points_2017 :
        print current_mass_point
        
f_2016.close()
f_2016 = open("./Zprime_to_NN_2016.txt", 'r')
print "=========== 2018 missed mass points ==========="
for line in f_2016:
    if not line: break
    line = line[0:-1]
    mass_point = line.split("/")[1]
    channel = mass_point.split("_")[1][0:2]
    if "Mu" in channel:
         channel = channel + channel
    mZp = mass_point.split("_")[2][6:]
    mN = mass_point.split("_")[3][1:]
    current_mass_point = channel + "_" + mZp + "_" + mN
    if current_mass_point not in all_mass_points_2018 :
        print current_mass_point
        
f_2016.close()
f_2017.close()
f_2018.close()

