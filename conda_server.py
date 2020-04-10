#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:18:07 2020

@author: Eric Chen, Graduate Center, CUNY

@ Prof. Kurtzman's Lab
"""

# Handling the conoserver PTM xml file

from Bio.SeqUtils import molecular_weight



import os

import numpy as np

import xlsxwriter


os.chdir("/Users/eric/Dropbox/Mande_projects/Tanya-project")

file=open("conoserver_protein.xml").readlines()

server_library=xlsxwriter.Workbook("Conno_server_library.xlsx")
library_sheet=server_library.add_worksheet()

library_sheet.write(0,0,"Record_ID")

library_sheet.write(0,1,"Record_seq")

server_mass=xlsxwriter.Workbook("Cono_server_MS.xlsx")

mass_sheet=server_mass.add_worksheet()

mass_sheet.write(0,0,"Record_ID")
mass_sheet.write(0,1,"Record_seq")
mass_sheet.write(0,2,"Max. MW")
mass_sheet.write(0,3,"PTM")





class cono:
    
    def __init__(self,seq,monomass,modification):
        
        self.seq=seq
        self.monomass=monomass
        self.modification=modification
        
cono_server=[]

    
output=open("conoserver.txt","w+")


i=0

s=0  # qualified sequence records from Conoserver

row=1 # put the record to excel

seq_lib=set()

for line in file:
    
    if "<sequenceModifications>" in line:
        #print(i)
        
        PTM_end=i + file[i:].index("</sequenceModifications>\n")
        
        seq=file[i-1].split(">")[1].split("<")[0]
        
        
        
        if ("." in file[PTM_end + 3].split(">")[1].split("<")[0]) and ("X"  not in seq) and ("x" not in seq):
            
            #mono_index=i + file[i:].index("*<monoisotopicMass>*")
            
            steps=[i for i , item in enumerate(file[i:]) if item.startswith("<monoisotopicMass>")]
            
            mono_index=i + steps[0]
            
            mono_M=float(file[mono_index].split(">")[1].split("<")[0])
        
            mod=file[i+1:PTM_end]
            
            sorted_seq=''.join(sorted(seq))
            
            
            t=cono(seq,mono_M,mod)
            
            #print("Seq: {0}, mono_mass: {1}, PTM: {2}".format(seq, mono_M, mod))
            
            cyc_num=seq.count("C")
            
            calculated_mass=molecular_weight(seq,"protein",monoisotopic=True) - 2* int(cyc_num/2)
            
           
        
        
        
            out="Seq: {0}, mono_mass: {1}, calculated_mass: {2}, diff: {3}, PTM: {4}".format(seq,mono_M,calculated_mass,calculated_mass-mono_M,mod)
            
            #only pull out th Carboxylic E, Hydro-proline, Bromide-W modification
            screen=[(("Gla"in i) or ("O" in i) or ("BTr" in i))  for i in mod]
            
            if np.all(np.array(screen)==1) and sorted_seq not in seq_lib:
                print(out+"\n")
            
            
            
                output.write("####seq{0}####\n".format(s))
                output.write(out+"\n")
                cono_server.append(t)
                
                library_sheet.write(row,0,s)
                library_sheet.write(row,1,seq)
                
                mass_sheet.write(row,0,s)
                mass_sheet.write(row,1,seq)
                mass_sheet.write(row,2,mono_M)
                
                simpled_mod=[i.split()[2].split("=")[-1] for i in mod]
                
                mass_sheet.write(row,3,str(simpled_mod))
                
                
                seq_lib.add(sorted_seq)
                
                s+=1
                row+=1
            
          
    
    i=i+1
        
print(len(cono_server))    

output.close() 

server_library.close()
server_mass.close() 

#molecular_weight("CCDDSECSTSCWPCCY","protein",monoisotopic=True)  

    
