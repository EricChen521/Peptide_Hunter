#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:04:53 2019

@author: Eric Chen, Graduate Center, CUNY

@ Prof. Kurtzman's Lab
"""

import os

#print(os.getcwd())

import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import argparse

import pandas as pd


parser=argparse.ArgumentParser(description="Welcome to peptide hunter!")

parser.add_argument('-f', type=str, help="The file name of peptide library, file formart should be xlxs " )

parser.add_argument('-P', type=int, help="The max number of phosphorylation allowed, 1000 represents all possible phosphorylation", default=0)

parser.add_argument("-S", type=int, help="The max number of disulfide bond allowed, 1000 represents all possible disulfide bond",default=0)

parser.add_argument("-t",type=float,help="The tolerance of Mass Spec accurancy,default=0.5 Dalton",default=0.002)
args=parser.parse_args()


# count the how many different PTM added

K=int(bool(args.P))+int(bool(args.S))

# read the peptide library xlsx file toxin_output.xlsx



library=pd.ExcelFile("all_toxin_no_duplicate.xlsx")

#print(library.sheet_names)

sheet_name=library.sheet_names[0]


df=library.parse(sheet_name)

#print(df["Record_seq"])

# df.shape (1151,7)

# read the record ID and Molecular weight to dictonary 
Database_1=[]


# Include the special modification 
Database_2=[]



from Bio.SeqUtils import molecular_weight

#define a peptide class

class Peptide:
    def __init__(self,ID,seq,MW):
        self.name=ID
        self.seq=seq
        self.MW=MW
        
        
for i in range(df.shape[0]):
    ID=df.iloc[i]["Record_ID"]
    seq=df.iloc[i]["Record_seq"]
    cys_num=seq.count("C")
    
    # calculate the original pepyide MW
    MW=float("%0.3f" % molecular_weight(seq,"protein"))
    

    Database_1.append(Peptide(ID,seq,MW))
    
    
    #N-terminal 
    
    #pyroglutamate Delta_Mass= -17 in GLn(Q), Delta_mass=-18 in Glu(E)
    if seq[0]=="Q": 
        MW=MW-17
    if seq[0]=="E":
        MW=MW-18
    
    # C-terminal 
    
    #  C-terminal Gly amidation Delta_Mass=-1
    
    if seq[-1]=="G":
        MW=MW-58
    
    # C-terminal CDPGLL Delata_Mass=-301.3818(GLL)-17(OH)+16(NH2)=-302.3818
    if seq[-6:]=="CDPGLL":
        MW=MW-302.3818
    
    Database_2.append(Peptide(ID,seq,MW))
    
    

 


# Calculate the maxium Phosphorylation sites of all peptides 

P_num=0  # Phosphorylation Site
S_num=0  # Cys-Cys disulfar bond

Pro_num=0 # hydroxyProline Delta_Mass=+16
Glu_num=0 #GamahydroxyGlu Delta_Mass=+16   ##how to distinguish these two 



original_MW_1=[]
 
  
for i in Database_1:
    
    # Serine(S), Threonine(T), Tyrosine(Y)
    P_num=i.seq.count("S")+i.seq.count("T")+i.seq.count("Y")
    
    # Cys 
    S_num=i.seq.count("C")
    
    # Pro
    
    Pro_num=i.seq.count("P")
    
    # Glu
    
    Glu_num=i.seq.count("E")
    
    
    original_MW_1.append(i.MW)
        

original_MW_2=[]

for i in Database_2:
    original_MW_2.append(i.MW)


# First dimension labels the peptide index in Database

Peptide_idx=df.shape[0]

import numpy as np

#args.P=3
#args.S=2

Dimension=(Peptide_idx,min(args.P,P_num)+1,min(args.S,S_num)+1,Pro_num+1,Glu_num+1)

#Dimension=(Peptide_idx,2,3)

MW_values=np.zeros(Dimension)

MW_0=np.array(original_MW_1) # do not in include spefical modification

MW_values[:Peptide_idx,0,0,0,0]=MW_0

#print(MW_values)

for d1 in range(1,Dimension[1]): # d1 dimension is phosphorylation
    MW_values[:Peptide_idx,d1,0,0,0]=MW_0 + 80*d1
  
MW_1=MW_values[:Peptide_idx,:Dimension[1],0,0,0]

for d2 in range(1,Dimension[2]):
    MW_values[:Peptide_idx,:Dimension[1],d2,0,0]=MW_1-2*d2

MW_2=MW_values[:Peptide_idx,:Dimension[1],:Dimension[2],0,0]

for d3 in range(1,Dimension[3]):
    MW_values[:Peptide_idx,:Dimension[1],:Dimension[2],d3,0]=MW_2+16*d3
    
MW_3=MW_values[:Peptide_idx,:Dimension[1],:Dimension[2],:Dimension[3],0]

for d4 in range(1,Dimension[4]):
    MW_values[:Peptide_idx,:Dimension[1],:Dimension[2],:Dimension[3],d4]=MW_3 + 16*d4
    


# inlcude special terminal modification
    
    
MW_values_01=np.zeros(Dimension)

MW_01=np.array(original_MW_2) # do not in include spefical modification

MW_values_01[:Peptide_idx,0,0,0,0]=MW_01

#print(MW_values)

for d1 in range(1,Dimension[1]): # d1 dimension is phosphorylation
    MW_values_01[:Peptide_idx,d1,0,0,0]=MW_01 + 80*d1
  
MW_1=MW_values_01[:Peptide_idx,:Dimension[1],0,0,0]

for d2 in range(1,Dimension[2]):
    MW_values_01[:Peptide_idx,:Dimension[1],d2,0,0]=MW_1-2*d2

MW_2=MW_values_01[:Peptide_idx,:Dimension[1],:Dimension[2],0,0]

for d3 in range(1,Dimension[3]):
    MW_values_01[:Peptide_idx,:Dimension[1],:Dimension[2],d3,0]=MW_2+16*d3
    
MW_3=MW_values_01[:Peptide_idx,:Dimension[1],:Dimension[2],:Dimension[3],0]

for d4 in range(1,Dimension[4]):
    MW_values_01[:Peptide_idx,:Dimension[1],:Dimension[2],:Dimension[3],d4]=MW_3 + 16*d4
        

# read the CESI_MS data

file1=pd.ExcelFile("LCMS.xlsx")



import sys

sys.stdout=open("output.txt","w+")
file1_sheets=file1.sheet_names

for species in file1_sheets:
    print("In the venom: "+species+", the detected peptides are:")
    
    sheet=file1.parse(species)
     
    #mw_values=[]
    for i in range(sheet.shape[0]):
        
        
        value=float("%0.3f" % sheet.iloc[i]["Max. MW"])
        
        tolerance=0.5
        tolerance=args.t
        coordinates_0=np.argwhere(abs(MW_values-value)<tolerance) 
        coordinates_1=np.argwhere(abs(MW_values_01-value)<tolerance)
        
        
        if len(coordinates_0):
            print("################### Do not include terminal modification ################")
            for position in coordinates_0:
                peptide_index=position[0]
                P_site=position[1]
                S_bond=position[2]
                Pro_hydroxy=position[3]
                Glu_hydroxy=position[4]
                peptide=Database_1[peptide_index].seq
                P_max=peptide.count("S")+peptide.count("T")+peptide.count("Y")
                Cys_max=peptide.count("C")
            
                if P_site<=P_max and S_bond <=int(0.5 * Cys_max):
                
                    print("Record_ID: {0}, Seq: {1}, Detected MW: {2}".format(Database_1[peptide_index].name,Database_1[peptide_index].seq,value))
                    print("{0} site(s) are phosphorylated".format(P_site))
                    print("{0} Cys-Cys disulfar bond(s) are formed".format(S_bond))
                    print("{0} Proline are hydroxygened".format(Pro_hydroxy))
                    print("{0} Glu are hydroxygened".format(Glu_hydroxy))
                
        if len(coordinates_1):
            
            print("################### Include terminal modification ################")
            
            for position in coordinates_1:
                peptide_index=position[0]
                P_site=position[1]
                S_bond=position[2]
                Pro_hydroxy=position[3]
                Glu_hydroxy=position[4]
                peptide=Database_1[peptide_index].seq
                P_max=peptide.count("S")+peptide.count("T")+peptide.count("Y")
                Cys_max=peptide.count("C")
                if P_site<=P_max and S_bond <=int(0.5 * Cys_max):
                
                    print("Record_ID: {0}, Seq: {1}, Detected MW: {2}".format(Database_1[peptide_index].name,Database_1[peptide_index].seq,value))
                    print("{0} site(s) are phosphorylated".format(P_site))
                    print("{0} Cys-Cys disulfar bond(s) are formed".format(S_bond))
                    print("{0} Proline are hydroxygened".format(Pro_hydroxy))
                    print("{0} Glu are hydroxygened".format(Glu_hydroxy))
