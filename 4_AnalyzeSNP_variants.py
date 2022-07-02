#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 16:51:14 2020

@author: minapeyton
"""
######## PART 1: Analysis of single-nucleotide variations in a population
# 1. Open and read in chromosomal sequences in slim_chr2_seq.fasta for 
# each of the 90 individuals. Store identifier and sequence information 
# in separate lists, one element in each list per individual.


# Process one line of identifier / sequence information

# Open the sequence file in read mode.
infile = open("slim_chr2_seq.fasta", 'r')

# Open two empty files in write mode
outfile = open("identifier.txt", 'wt')
print("Identifier",end='\n', file = outfile)

outfile2 = open("seqinfo.txt", 'wt')
print("Sequence Information",end='\n', file = outfile2)

##### Code for processing the first line only
# Read in the first line of the file
line1 = infile.readline()

# Remove any trailing newline characters from the sequence
line1 = line1.rstrip()

# find > identifier and store in namelist
line1.find(">")
name = line1.split()
namelist = []
namelist.append(name[1])

# find sequence information and store is seqlist
line1 = infile.readline()
seq =line1.split()
seqlist = []
seqlist.append(seq[0])

# print to outfile
print(namelist,sep='\,', end='\n', file = outfile)
print(seqlist,sep='\,', end='\n', file = outfile2)

print("Done")

# Close the input file -- (NOTE this is not inside the for loop, so it should not be indented)
infile.close()
# Close the output file
outfile.close()

#############################################
##############################################

# Process whole "slim_chr2_seq.fasta" file

# Open the sequence file in read mode.
infile = open("slim_chr2_seq.fasta", 'r')

# Open an empty file called 'aaoutput.txt' in write mode
outfile = open("identifier.txt", 'wt')
print("Identifier",end='\n', file = outfile)

outfile2 = open("seqinfo.txt", 'wt')
print("Sequence Information",end='\n', file = outfile2)

# create empty list to store identifier & sequence info

namelist = []
seqlist = []

##### for loop to process whole file and create two lists for each individual

for line in infile:
    if line.find(">") == 0:
        line = line.rstrip()
        name = line.split()
        namelist.append(name[1])
    else:
        seq =line.split()
        seqlist.append(seq[0])

# print to outfile
print(namelist,sep='\,', end='\n', file = outfile)
print(seqlist,sep='\,', end='\n', file = outfile2)

print("Done")

# Close the input file -- (NOTE this is not inside the for loop, so it should not be indented)
infile.close()
# Close the output file
outfile.close()
outfile2.close()

################################################

# 2. Open and read in SNP information stored in slim_chr2_SNPS.vcf 
# (which can be opened as a text file). Store the information in this file 
# in seven different lists, one for each column of data (open the file in a 
# text editor before writing your code to understand its format). Remember 
# to close any files you open!

# Open the sequence file in read mode.
infile2 = open("slim_chr2_SNPS.vcf", 'r')

# Open seven empty files in write mode
outfile3 = open("chr_list.txt", 'wt')
# print("Chromosome",end='\n', file = outfile3)

outfile4 = open("rsid_list.txt", 'wt')
# print("rsid",end='\n', file = outfile4)

outfile5 = open("slim_pos_list.txt", 'wt')
# print("slim_pos",end='\n', file = outfile5)

outfile6 = open("genome_pos_list.txt", 'wt')
# print("genome_pos",end='\n', file = outfile6)

outfile7 = open("ref_list.txt", 'wt')
# print("ref",end='\n', file = outfile7)

outfile8 = open("alt_list.txt", 'wt')
# print("alt",end='\n', file = outfile8)

outfile9 = open("geneinfo_list.txt", 'wt')
# print("geneinfo",end='\n', file = outfile9)

# open seven empty lists to read and store SNP information

Chr_list = []
rsid_list = []	
slim_pos_list = []	
genome_pos_list =[]	
ref_list = []	
alt_list = []
geneinfo_list = []

##### For loop to process file

for line in infile2:
    line = line.rstrip()
    line = line.split(sep='\t')
    Chr_list.append(line[0])
    rsid_list.append(line[1])
    slim_pos_list.append(line[2])
    genome_pos_list.append(line[3])
    ref_list.append(line[4])
    alt_list.append(line[5])
    geneinfo_list.append(line[6:])

# print to outfile
print(Chr_list,sep='\,', end='\n', file = outfile3)
print(rsid_list,sep='\,', end='\n', file = outfile4)
print(slim_pos_list,sep='\,', end='\n', file = outfile5)
print(genome_pos_list,sep='\,', end='\n', file = outfile6)
print(ref_list,sep='\,', end='\n', file = outfile7)
print(alt_list,sep='\,', end='\n', file = outfile8)
print(geneinfo_list,sep='\,', end='\n', file = outfile9)

print("Done")

# Close the input file -- (NOTE this is not inside the for loop, so it should not be indented)
infile2.close()
# Close the output file
outfile3.close()
outfile4.close()
outfile5.close()
outfile6.close()
outfile7.close()
outfile8.close()
outfile9.close()

#########################################
# 3. For each provided variant, calculate the frequency of that SNP in 
# the population of 90 individuals. In other words, how many of the 90 
# individuals have the corresponding base at the specified location for 
# that variant? Store these variant frequencies in a separate list.

pos = slim_pos_list[1:]
ref = ref_list[1:]
alt = alt_list[1:]
seqlist_copy = seqlist[1:]

# create SNP alternative dictionary (position, alt_list)
alt_dict = {}

for i in range(len(pos)):
    pos2 = pos[i]
    alt2 = alt[i]
    alt_dict[int(pos2)+1] = alt2
    
# create ref dictionary (position, ref_list)

ref_dict = {}

for i in range(len(ref)):
    pos2 = pos[i]
    ref2 = ref[i]
    ref_dict[int(pos2)+1] = ref2    

# create new list to store variant frequencies in a list 
snp_frequencies = {}

ref_frequencies = {}

# loop through keys of alt.dictionaries to count the frequency of that SNP 
# in the 90 individuals

for key in alt_dict:
    for seq in seqlist_copy:
        c = list(seq)
        if c[key] == alt_dict[key] and key in snp_frequencies:
            snp_frequencies[key] = snp_frequencies[key] + 1
        elif c[key] == alt_dict[key]:
                snp_frequencies[key] = 1
        else:
            continue
            
            
for (key,val) in snp_frequencies.items():
		print(key,"-", val)
         
for key in ref_dict:
    for seq in seqlist_copy:
        c = list(seq)
        if c[key] == ref_dict[key] and key in ref_frequencies:
            ref_frequencies[key] = ref_frequencies[key] + 1
        elif c[key] == ref_dict[key]:
                ref_frequencies[key] = 1
        else:
            continue

for (key,val) in ref_frequencies.items():
		print(key,"-", val)
        

for seq in seqlist_copy:
    c = list(seq)
    for key in ref_dict:
        if c[key] == ref_dict[key] and key in ref_frequencies:
            ref_frequencies[key] = ref_frequencies[key] + 1
        elif c[key] == ref_dict[key]:
                ref_frequencies[key] = 1
        else:
            continue
    print()

for (key,val) in ref_frequencies.items():
		print(key,"-", val)
        
for seq in seqlist_copy:
    c = list(seq)
    for key in alt_dict:
        if c[key] == alt_dict[key] and key in snp_frequencies:
            snp_frequencies[key] = snp_frequencies[key] + 1
        elif c[key] == alt_dict[key]:
                snp_frequencies[key] = 1
        else:
            continue
    print()

for (key,val) in snp_frequencies.items():
		print(key,"-", val)

#########################################
# 4. Find and print the ids of the SNPs in the population that occur 
# with the highest and lowest frequency.
  




#########################################
# 5. Open a new output file called slim_chr2_SNPs_withfrequencies.vcf and 
# print the following information for each SNP:
# a. Chromosome Number
# b. SNP Name
# c. SNP “Slim” Position
# d. SNP Reference SNP
# e. SNP Alternative SNP
# f. SNP Frequency
# g. Gene Disease Information  
    

outfile10 = open("slim_chr_SNPs_withfrequencies.vcf", 'wt')
# print("Chr_list", "SNP_name", "slim_pos_list", ref_list", "alt_list", "snp_frquencies"end='\n', file = outfile9)


##########  Part II: Characterizing SNPs
# 
    
        
















