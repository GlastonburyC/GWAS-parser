#!/bin/python
import re
import os
import sys

try:
    gwas = sys.argv[1]
except:
    print "Usage: python extract.py [GWAS Catalog]"


handle = file(gwas,"r")
handle_result = file("result.txt","w")
line = handle.readlines()
for entry in line[1:-1]: 
    entry_split = entry.split('\t')
    p = entry_split[27]
    
    try:
        p = float(p)
        if p<5.00E-8:
            m = re.search('(?<=-)\w',entry_split[20])
            allele = m.group()
            if allele == "?":
                allele = [NR]
           
            handle_result.write(entry_split[13]+'\t'+entry_split[7]+'\t'+entry_split[27]+'\t'+entry_split[30]+ '\t'+ entry_split[31] +'\t'+entry_split[21]+'\t'+allele+'\t'+entry_split[26]+'\t'+entry_split[1]+'\n')
    except:
        pass
        
handle.close()

handle = file("result.txt","r")
dict_trait = {}
content = handle.readlines()
for line in content:
    entry = line.rstrip('\n').split('\t')
    if not(entry[1] in dict_trait):
        dict_trait[entry[1]] = []
    #Remove 4 underscores from OR/Beta value if it exists
   #Can be removed if NHGRI catalog fixes the spacing
    if entry[3][0]==" ":
        modentry3 = entry[3].strip('    ')
        dict_trait[entry[1]].append( (entry[0],entry[1],entry[2],modentry3,entry[4],entry[5],entry[6],entry[7],entry[8])  )
    else:
        dict_trait[entry[1]].append( (entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],entry[6],entry[7],entry[8])  )

os.remove("result.txt")



handle_index = file("index","w")
counter = 0
path ='./traitFiles/'
if not os.path.exists(path):
    os.makedirs(path)
header = 'Gene' +'\t' + 'Trait' +'\t' + 'p_Value' +'\t' + 'OR/Beta'+ '\t' '95%CI'+'\t' + 'SNPrsID' +'\t' + 'Risk_Allele' +'\t' + 'Risk_Allele_Freq' +'\t' +'PubMedID' +'\n'
for key in dict_trait.keys():
    handle = file(path+str(counter),"w")
    templist = dict_trait[key]
    handle.write(header)
    for tuple in templist:
        for string in tuple:
            new_string = string.replace (" ","_")
            handle.write(str(new_string)+'\t')
            
        handle.write('\n')
    handle.close()
    handle_index.write(str(counter)+'\t'+key+'\n')
    counter += 1
   


handle_index.close()
