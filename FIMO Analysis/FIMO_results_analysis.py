# -*- coding: utf-8 -*-
"""
Created on Tue May 18 09:15:12 2021

@author: squir
"""

#setwd C:\Users\squir\OneDrive\Documents

#open the main datafile
FIMO_Dofs = open("fimo_Dofs.tsv", 'r')
FIMO_TCPs = open("fimo_TCP_Thresh.001.tsv", 'r')
FIMO_TCPs_strict = open("fimo_TCP.tsv", 'r')
FIMO_HDZIPs = open("fimo_HDzip.tsv", 'r')


#open the new file to deposit MEME formatted stuff
outfile = open("Hist_Dofs.csv", 'w')

#Create Header for MEME file
print("Nucleotide", "Dof_Hits", "TCP_Hits", "HD_ZIP_Hits","TCP_strict_Hits", sep = ",", end = "\n", file = outfile)

#Start a dictionary for each TF
Dict_dofs = {}
Dict_TCPs = {}
Dict_hdzip = {}
Dict_TCP_strict = {}
#create a key for each n.tide position, start counters at 0
for i in range(1,1025):
    Dict_dofs[i] = 0
    Dict_TCPs[i] = 0
    Dict_hdzip[i] = 0
    Dict_TCP_strict[i] = 0
    
#Loop through each line of the input file
for line in FIMO_Dofs:
    #rstrip each line
    line = line.rstrip()
    #skip empty lines and header lines
    if line != '' and line[0] == 'M':
        #split line on tabs
        temp = line.split("\t")
        #Add up hits at each positions
        for i in range(int(temp[3]),int(temp[4])+1):
            Dict_dofs[i] += 1

#Loop through each line of the input file
for line in FIMO_TCPs:
    #rstrip each line
    line = line.rstrip()
    #skip empty lines and header lines
    if line != '' and line[0] == 'M':
        #split line on tabs
        temp = line.split("\t")
        #Add up hits at each positions
        for i in range(int(temp[3]),int(temp[4])+1):
            Dict_TCPs[i] += 1

#Loop through each line of the input file
for line in FIMO_TCPs_strict:
    #rstrip each line
    line = line.rstrip()
    #skip empty lines and header lines
    if line != '' and line[0] == 'M':
        #split line on tabs
        temp = line.split("\t")
        #Add up hits at each positions
        for i in range(int(temp[3]),int(temp[4])+1):
            Dict_TCP_strict[i] += 1

#Loop through each line of the input file
for line in FIMO_HDZIPs:
    #rstrip each line
    line = line.rstrip()
    #skip empty lines and header lines
    if line != '' and line[0] == 'M':
        #split line on tabs
        temp = line.split("\t")
        #Add up hits at each positions
        for i in range(int(temp[3]),int(temp[4])+1):
            Dict_hdzip[i] += 1

for key in Dict_dofs:
    print(key, Dict_dofs[key], Dict_TCPs[key], Dict_hdzip[key], Dict_TCP_strict[key], sep = ",", end = "\n", file = outfile)

FIMO_Dofs.close()
FIMO_TCPs.close()
FIMO_HDZIPs.close()
FIMO_TCPs_strict.close()
outfile.close()
