#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:21:31 2019

@author: Eric Chen, Graduate Center, CUNY

@ Prof. Kurtzman's Lab
"""

import os 


import pandas as pd

import numpy as np

from Bio.SeqUtils import molecular_weight

import math

import matplotlib.pyplot as plt


os.chdir("/Users/eric/Desktop/Mande")



"""

This script only consider Carbamidomethyl( Fixed) and Amidated(variable) and Oxidation (Variable) 

and the signal with charges +1, +2 , +3 to match the library 


"""

#Step 1, transform the fasta txt file into excel file

file_txt=open("All_Putative_Toxins_WoD.fasta").readlines()



n_record=0

data=[]

for line in file_txt:
    
    
    
    if ">" in line:
        ID=line.split(">")[1]
        
        data.append([ID," "]) # induce a empty string here
        
        #worksheet.write(int(n_record/1),0,ID)
        
        
        
        n_record+=1
        
        
    else:
        
        
        seq=line.strip()
        
        data[n_record-1][1]+=seq
        
        
        #worksheet.write(int(n_record/1),1,seq)
    
seqs=[]

no_duplicate_data=[]

for i in range(len(data)):
    
    
    
    data[i][1]=data[i][1].strip()  # get rid of the empty string at the head of the seqence
    
    if data[i][1] not in seqs:
        
        seqs.append(data[i][1])
        
        no_duplicate_data.append(data[i])
        
  


df=pd.DataFrame(no_duplicate_data,columns=["Record_ID","Record_seq"])

df.to_excel("test_2.xlsx",index=False)


library=[]

class Peptide:
    
    
    
    def __init__(self,ID,seq,MW):
        self.name=ID
        self.seq=seq
        self.MW=MW
        #self.charge=charge
 



dir(df)

df.columns


for i in range(df.shape[0]):
    
    ID=df.iloc[i]["Record_ID"].strip()
    
    seq=df.iloc[i]["Record_seq"]
    
    #MW=float("%0.3f" % molecular_weight(seq,"protein",monoisotopic=True))
    MW=float("%0.3f" % molecular_weight(seq,"protein"))
    
    
    
    #charge=df.iloc[i]["Max. z"]
    
    cys_num=seq.count("C")
    
    
    ## fixed carbamidomethal modification on all cys 
    
    # + 57.021464 Da reference: parameters from mascot
    
    
    #MW= MW + (57.021464 * cys_num)
    
    library.append(Peptide(ID,seq,MW))

#print(library)

## Variable PTM modications, we only only one amidation and one oxidation 
    
# 1) Amidation -0.984016
    
# 2) Oxidation [ + 15.994915] Refernce https://www.ncbi.nlm.nih.gov/pubmed/8374164
    
# 3) disulfide # limit up to three
    
    
Peptide_idx=df.shape[0]

#print(Peptide_idx)

s_bonds=5

Dimension=(Peptide_idx,2,2,10)

MW_values=np.zeros(Dimension)

original_MW=[i.MW for i in library]

MW_values[:Peptide_idx,0,0,0]=np.array(original_MW)

MW_values[:Peptide_idx,1,0,0]=MW_values[:,0,0,0]- 0.984016 ## single amidation

MW_1=MW_values[:,:,0,0]

MW_values[:Peptide_idx,:,1,0]=MW_1 + 15.994915 ## single oxidation

# disulfer bonds

Mw_2=MW_values[:,:,:,0]

for s in range(1,s_bonds):
    
    MW_values[:,:,:,s]=Mw_2 - 2 *s



pure_m1=2695.9442
hit1_index=np.argwhere(np.absolute((MW_values-pure_m1))<0.1)

hit1_mass=molecular_weight("FDYESLWDTV","protein",monoisotopic=True)

print(hit1_mass)

print(hit1_index)


hit1=library[hit1_index[0][0]]

hit1.MW





pure_m2=1431.4793
hit2_index=np.argwhere(np.absolute((MW_values-pure_m2))<0.1)

print(hit2_index)

hit2=library[hit2_index[0][0]]

hit2.name


MW_values[324,0,0,3]

MW_values[324,0,0,0]

MW_values[324,1,0,2]

#MW_2=MW_values[:,:,:,0]
"""
for d_cys in range(1,10):
    
    MW_values[:,:,:,d_cys]=MW_2 - 2*d_cys

## Possible Mass value were filled 

## read the Mass spec data
"""

#analysis the mascot data 

file1=pd.ExcelFile("Peptide_Hunter/MascotResults.xlsx")

mascot_match={}

mascot_mass=[]

mascot_record=set()


file1_sheets=file1.sheet_names

for sheet in file1_sheets:
    
    contents=file1.parse(sheet)
    
    peptide_set=set()
    
    species=sheet
    
    for i in range(58,contents.shape[0]):
        
        record=contents.iloc[i][1]
        
        mass=contents.iloc[i][14]
        
        peptide_set.add(record)
        
        mascot_record.add(record)
        
        mascot_mass.append(mass)
        
    
    mascot_match.update({species:peptide_set})

print(mascot_record)

print(len(mascot_record))

print(mascot_mass)
    











#print("Haha")

LC_mass=[]
#CESI_mass=[]


LC_record=set()

LC_charge=[]


file2=pd.ExcelFile("CESI-MS.xlsx")

#file2=pd.ExcelFile("Pure_peptide.xlsx")

#file2=pd.ExcelFile("Peptide_Hunter/LCMS.xlsx")

file2_sheets=file2.sheet_names

for sheet in file2_sheets:
    
    contents=file2.parse(sheet)
    
    for i in range(contents.shape[0]):
        
        
        value=float(contents.iloc[i]["Max. MW"])
        
        charge=contents.iloc[i]["Max. z"]
        
        if not math.isnan(value):
                        
            LC_charge.append(int(charge[0]))
            
            
            LC_mass.append(value)
        
print(LC_mass)

#show the distribution of mass detected by tamdem mass and lC_MS
        
plt.hist(mascot_mass,bins=10,normed=True,color="orange",alpha=0.5,label="Tandem MS")

plt.hist(LC_mass,bins=10,normed=True, color="blue",alpha=0.5,label="CESI_MS")

plt.hist(MW_values[:,0,0,0],normed=True, bins=10,alpha=0.5,color="grey",label="Transcriptome library")


plt.xlabel("peptide mass")
plt.ylabel("Normaled Frequency")

plt.legend()

plt.show()






#match=open("CESI_MS_signal_monoisopotic.txt","w+")

match=open("CESI_MS_match.txt","w+")

#match=open("LC_MS_Match_monoisotopic.txt","w+")

#match=open("LC_MS_match.txt","w+")

# set the peptide tolerance 0.1

tolerance=0.1

LC_MS_match_mass=[]

LC_MS_match_charge=[]

LC_MS_match_record=[]

for species in file2_sheets:
    
    print("In the venom: " + species+", the detected peptides are:",file=match)
    
    sheet=file2.parse(species)
        
    for i in range(0,sheet.shape[0]):
        
        
        
        value=float(sheet.iloc[i]["Max. MW"])
        
        if not math.isnan(value):
            
            
        
            charge=sheet.iloc[i]["Max. z"]
        

            coordinates_1=np.argwhere(np.absolute((MW_values-value))<tolerance)
            
        
            if len(coordinates_1):
                
                LC_MS_match_mass.append(value)
                
                LC_MS_match_charge.append(int(charge[0]))
                
                
                
        
                for position in coordinates_1:
                    #print(position)
            
                    Peptide_index=position[0]
                    Amidation_num=position[1]
                    Oxidation_num=position[2]
                
                    S_num=position[3]
            
                    peptide=library[Peptide_index]
                    peptide_seq=peptide.seq
                    
                    
                    
                    peptide_cys_num=peptide.seq.count("C")
                    
               
                
                    if  peptide_cys_num >= 2 * S_num:
                        
                        LC_MS_match_record.append(peptide.name)
                    
            
                        print("Record_ID: {0}, Seq: {1}, Detected MW: {2}, charges: {3}. Coordinate : {4}".format(peptide.name,peptide.seq,
                      value,charge,position),file=match)
    
                
match.close()



# show the correlation of mass with charge

size1=[5 for i in LC_mass]
plt.scatter(LC_mass, LC_charge,c="grey",s=size1,alpha=0.5, label="All signal")

size2=[5 for i in LC_MS_match_mass]

plt.scatter(LC_MS_match_mass, LC_MS_match_charge,s=size2, alpha=0.5, color="red", label="Matched hits")

plt.xlabel("Peptide mass (Da)")

plt.ylabel("Peptide charge (+)")

plt.legend(loc="center right")

plt.show()


## compaire the hit from LC_MS and tandem MS

print(LC_MS_match_record)

print(len(LC_MS_match_record))

print(mascot_record)

print(len(set(mascot_record)))

print(len(set(LC_MS_match_record)))

from matplotlib_venn import venn2, venn2_circles

venn2([set(LC_MS_match_record),set(mascot_record)],set_labels=("CE_MS","Tandem MS"))











hooked_pep_set=set(hooked_peptide)

len(hooked_pep_set)



hooked_pep_ls=list(hooked_pep_set)


for i in library:
    if i.name in hooked_pep_set:
        
        print(i.MW)
        
        Tandom_MS_values.append(i.MW)
    

print(Tandom_MS_values) 

print(LC_MS_values)  



print(file1_sheets)

print(file2_sheets)






import matplotlib.pyplot as plt




plt.hist(Tandom_MS_values,bins=10,color="blue", label="Tandem MS")

plt.hist(LC_MS_values,bins=10,color="red", label="LC MS")

plt.legend()

plt.show()

Transcriptome_mass=[i for i in MW_values[:,0,0] ]  



bins=np.linspace(0,16000,50)
plt.hist(Transcriptome_mass,bins,alpha=0.7,label="Transcriptome Mass")
plt.hist(LC_mass,bins,alpha=0.5,label="LC-MS Mass")
plt.hist(CESI_mass,bins,alpha=0.7,label="CESI Mass")

plt.xlabel("Mass")
plt.ylabel("Counts")



plt.legend()

plt.savefig("Mass_distrubution_transcriptome_detected.png",dpi=400)

plt.show()


#print("hha")

 

     
#print(MW_values[:20,0,0]) 

#library              
          
#print(MW_values[:5,0,0])

#print(MW_values[:5,1,0])

#print(MW_values[:5,1,1])

                                 
                          
        
        
    
    
    
    
    

    
    

    
    
    
    



    

