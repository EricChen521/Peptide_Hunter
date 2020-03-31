#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:04:53 2019

@author: Eric Chen, Graduate Center, CUNY

@ Prof. Kurtzman's Lab
"""
# Add the post modificatio based on https://www.ncbi.nlm.nih.gov/pubmed/16314929


"""
In this script we will consider 

Fixed modification:
    
 1) Disulfide Bridge between Cys-Cys: delta_mass = -2 

Variable modification:

 2) Phosphorylation (S/T/Y): Delta_mass = +80
    
 3) Hydroxylation of Proline( P), Lysine(Lys,K): delta_mass = +16
    
 4) Amidation of C-terminal -XG(-75.03+1+16),  -XGK(-203.13+17), -XGR(-231.13+17), -XGRR(-387.23+17) to  -X-NH2: delta_mass= -58.03, -186.13,-214.13, -370.23
    Gly(75), R(174.2), K(146.2), based on https://web.expasy.org/peptide_mass/
    
 5) Carboxylation of glutamic acid (E), delata_mass = +44
    
 6) Bromination of Try(W): detlta_mass = +80-1= + 79
    
 7) Cyclization of N-terminal glutamine(Q), delta_mass = -17
  
"""


import sys

import numpy as np

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="Welcome to peptide hunter!")

parser.add_argument("-file",help="Mass Spec xlsx result file")
parser.add_argument("-P",type=int,help="Max number of phosphorylation allowed on S/T/Y, Default = 0;",default=0)
parser.add_argument("-OH",type=int,help="Max num of hydroxylation allowed on P/K, Default = 0;",default=0)
parser.add_argument("-CO2",type=int,help="Max num carboxylation allowed on E, Default = 0;",default=0)
parser.add_argument("-Br",type=int,help="Max num bromination allowed on W , Default = 0;",default=0 )
parser.add_argument("-NH2",help="Allow amidation on -XG,-XGR.-XGK,-XGRR if specified ",action="store_true")
parser.add_argument("-Cyc",help="Allow cyclization of N-terminal Q if specified;",action="store_true")
parser.add_argument("-Iso",help="Allow isotopic +1 +2 +3 shift if specified",action="store_true")
parser.add_argument("-t",type=float,help="The tolerance of Mass Spec accurancy,Default=0.1 Dalton;",default=0.1)

args=parser.parse_args()


class PTM:
    def __init__(self,name,delta_mass,AA_sites):
        
        self.name=name
        self.delta_mass=delta_mass
        self.AA_sites=AA_sites

# count the how many different PTM added

#Mass_shift={"P":80,"OH":16,"NH2":[-58.03,-186.13,-214.13,-370.23],"CO2":44,"Br":79,"Cyc":-17}
"""
Selected_PTM=[]

if args.P:
    Selected_PTM.append(PTM("P",80,["S","T","Y"]))
if args.OH:
    Selected_PTM.append(PTM("OH",16,["P","K"]))
if args.NH2:
    Selected_PTM.append(PTM("NH2",[-58.03,-186.13,-214.13,-370.23],["G","GK,""GR","GRR"]))
if args.CO2:
    Selected_PTM.append(PTM("CO2",44,["E"]))
if args.Br:
    Selected_PTM.append(PTM("Br",79,["W"]))
if args.Cyc:
    Selected_PTM.append(PTM("Cyc",-17,["Q"]))
if args.Iso:
    Selected_PTM.append(PTM("Iso",[1,2,3],[""]))

"""


# read the peptide library xlsx file toxin_output.xlsx



library=pd.ExcelFile("toxins_output.xlsx")


df=library.parse("toxins_output")


# df.shape (1151,7)

# read the record ID and Molecular weight to dictonary 
Database_raw=[]



from Bio.SeqUtils import molecular_weight

molecular_weight("G","protein",monoisotopic=True)

#define a peptide class

class Peptide:
    def __init__(self,ID,seq,MW):
        self.name=ID
        self.seq=seq
        self.MW=MW
        
MW_raw=[] # raw pepyide mass without any modification       
for i in range(df.shape[0]):
    ID=df.iloc[i]["Record_ID"]
    seq=df.iloc[i]["Record_seq"]
    
    # calculate the original pepyide MW uisng monoisotopic mass
    
    MW=float("%0.3f" % molecular_weight(seq,seq_type="protein",monoisotopic=True))
    
    # take the fixed disulfide bridge 
    
    Cys_num=seq.count("C")
    
    MW= MW- 2* int(Cys_num/2)
    MW_raw.append(MW)
    

    Database_raw.append(Peptide(ID,seq,MW))
    


# First dimension labels the peptide index in Database

d1=df.shape[0]



#First dimention is the raw peptide
MW_1=np.zeros(d1)

MW_1=np.array(MW_raw)

#Second dimension is phosphorylation

d2=min(args.P,10) +1



MW_2=np.zeros((d1,d2))

MW_2[:d1,0]=MW_1

for i in range(1,d2):
    MW_2[:d1,i]=MW_1 + 80*i

# Third dimension is Hydroxylation 
    
d3=min(args.OH, 10) +1

MW_3=np.zeros((d1,d2,d3))

MW_3[:d1,:d2,0]=MW_2

for i in range(1,d3):
    
    MW_3[:d1,:d2,i]=MW_2+ 16*i
    
# Fourth dimension is Carboxylation 
    
d4=min(args.CO2,10) + 1

MW_4=np.zeros((d1,d2,d3,d4))

MW_4[:d1,:d2,:d3,0]=MW_3

for i in range(1,d4):
    MW_4[:d1,:d2,:d3,i] = MW_3 + i * 44

# Fifth dimension is Bromination 
    
d5=min(args.Br,10) +1 

MW_5=np.zeros((d1,d2,d3,d4,d5))

MW_5[:d1,:d2,:d3,:d4,0]=MW_4

for i in range(1,d5):
    
    MW_5[:d1,:d2,:d3,:d4,i] = MW_4 + i * 79

# Sixth dimension for Amidation in C-terminal
    

  
if args.NH2:
    
    d6=5

    MW_6=np.zeros((d1,d2,d3,d4,d5,d6))
    
    MW_6[:d1,:d2,:d3,:d4,:d5,0]=MW_5


    NH2_shift=[-58.03,-186.13,-214.13,-370.23]
    
    for i in range(1,d6):
    
        MW_6[:d1,:d2,:d3,:d4,:d5,i]=MW_5 + NH2_shift[i-1]

    
else:
    
    d6=1
    
    MW_6=np.zeros((d1,d2,d3,d4,d5,d6))
    

    
    MW_6[:d1,:d2,:d3,:d4,:d5,0]=MW_5




    

    
# Seventh dimensition for Cyclization 
    
if args.Cyc:
    
    d7=2
    
    MW_7=np.zeros((d1,d2,d3,d4,d5,d6,d7))
    
    MW_7[:d1,:d2,:d3,:d4,:d5,:d6,0]=MW_6
    
    MW_7[:d1,:d2,:d3,:d4,:d5,:d6,1]=MW_6 -17 
    
else:
    
    d7=1
    
    MW_7=np.zeros((d1,d2,d3,d4,d5,d6,d7))
    
    MW_7[:d1,:d2,:d3,:d4,:d5,:d6,0]=MW_6
    
    
        
# Eighth Dimention for Isotopic (+1 +2,+3)
    
if args.Iso:
    
    d8 = 4
    
    MW_8=np.zeros((d1,d2,d3,d4,d5,d6,d7,d8))
    
    MW_8[:d1,:d2,:d3,:d4,:d5,:d6,:d7,0]=MW_7
    
    for i in range(1,d8):
        
        MW_8[:d1,:d2,:d3,:d4,:d5,:d6,:d7,i]=MW_7 +i

else:
    
    d8=1
    
    MW_8=np.zeros((d1,d2,d3,d4,d5,d6,d7,d8))
    
    MW_8[:d1,:d2,:d3,:d4,:d5,:d6,:d7,0]=MW_7
    
    


# read the Mass Spec data

#file1=pd.ExcelFile("LCMS.xlsx")
    
#file1=pd.ExcelFile("Pure_peptides.xlsx")
    
file1=pd.ExcelFile(args.file)



sys.stdout=open("P["+str(args.P)+"]_"+"OH["+str(args.OH)+"]_CO2["+str(args.CO2)+"]_Br["+str(args.Br)+"]_NH2["+str(args.NH2)+"]_Cyc["+str(args.Cyc)+"]_Iso["+str(args.Iso)+"]_tol["+str(args.t)+"]_output.txt","w+")

file1_sheets=file1.sheet_names

for species in file1_sheets:
    print("In the venom: "+species+", the detected peptides are:")
    
    sheet=file1.parse(species)
     
    #mw_values=[]
    for i in range(sheet.shape[0]):
        
        value=float("%0.3f" % sheet.iloc[i]["Max. MW"])
        
        #tolerance=args.t
        tolerance=args.t
        
        coordinates=np.argwhere(abs(MW_8-value)<tolerance) 
        
        if len(coordinates):
            
            for position in coordinates:
                
                peptide_index=position[0]
                peptide_seq=Database_raw[peptide_index].seq
                
                P_max=peptide_seq.count("S") + peptide_seq.count("T") + peptide_seq.count("Y")
                P_num=position[1]
                
                OH_max=peptide_seq.count("P") + peptide_seq.count("K")  
                OH_num=position[2]
                
                CO2_max=peptide_seq.count("E")               
                CO2_num=position[3]
                
                Br_max=peptide_seq.count("W")                
                Br_num=position[4]
                
                
                NH2_site=position[5]
                
                Cyc_able=int(peptide_seq[0]=="Q")
                Cyc_site=position[6]
                
                Iso_position=position[7]
                
                
                if (P_num <= P_max) and (OH_num <= OH_max) and (CO2_num <= CO2_max) and (Cyc_site <= Cyc_able):
                    
                    if ((NH2_site == 0) or (NH2_site == 1 and peptide_seq[-1] == "G") or (NH2_site ==2 and peptide_seq[-2:]=="GK") or (NH2_site ==3 and peptide_seq[-2:]=="GR") or (NH2_site == 4 and peptide_seq[-3:]=="GRR")):

                    
                        print("Record_ID: {0}, Seq: {1}, Detected MW: {2}".format(Database_raw[peptide_index].name,Database_raw[peptide_index].seq,value))
                        
                        print("P_sites:{0}, OH_sites:{1}, CO2_sites:{2}, Br_site:{3}, NH2_cleavage: {4}, Cyc_sites:{5},Iso_plus:{6}".format(P_num,OH_num,CO2_num,Br_num,["None","G","GK,""GR","GRR"][NH2_site],Cyc_site,Iso_position))
                    
              
