# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 09:35:35 2021

@author: squir
"""
#wdir = 'N:/Byron', open raw data file
from tkinter import filedialog as fd
filename = fd.askopenfilename()
Seg_data = open(filename, 'r')
#Seg_data = open("Intensity_data_ATML1_crosses_v2.csv", 'r')

#open the new file to deposit csv formatted data
outfile = open(filename[:-4] + "_parsed" + ".csv", 'w')

#Create Header for parsed csv file
print("Image_name", "Genotype", "Plant_number","Cell_type", "Count", "Total_area", "Avg_area","Percent_area","Mean_intensity","Integrated_density", sep = ",", end = "\n", file = outfile)

for line in Seg_data:
    #rstrip each line
    line = line.rstrip()
    #skip empty lines
    if line == '':
        continue
    if line[0:4] == 'Name':
        continue
    #split line on tabs
    temp = line.split("\t")
    #save each original column
    Image_name = temp[0]
    Count = temp[1]
    Total_area = temp[2]
    Avg_area = temp[3]
    Percent_area = temp[4]
    Mean_intensity = temp[5]
    Integrated_density = temp[6]
    #Further split the file name and extract information
    temp2 = temp[0].split("_")
    Cell_type = temp2[1][0:5]
    temp3 = temp2[2].split("-")
    Genotype = temp3[0]
    Plant_number = temp3[2][3:]
    #Print information to new file
    print(Image_name, Genotype, Plant_number, Cell_type, Count, Total_area, Avg_area,Percent_area,Mean_intensity,Integrated_density, sep = ",", end = "\n", file = outfile)

Seg_data.close()
outfile.close()
