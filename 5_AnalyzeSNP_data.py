#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:15:17 2020

@author: minapeyton
"""

        
# Part 1: define function to check if string represents a valid Homo sapiens 
# genome locus identifier
        
def is_valid_humanlocus(string):
    ''' Validates whether string contains Homo sapien genome locus identifiers '''
    if '.' in string:
        if 'p' in string or 'q' in string:
            fragments = string.split('.')
            left_frag = fragments[0]
            subband = fragments[1]
            subband = int(subband)
            if 'p' in left_frag:
                frag2 = left_frag.split('p')
                chromosome = frag2[0]
                band = frag2[1]
                assert chromosome.isdigit()
                assert band.isdigit()
                #print (chromosome)
                #print (band)
                #print (subband)
            elif 'q' in left_frag:
                frag2  = left_frag.split('q')
                chromosome = frag2[0]
                band = frag2[1]
                assert chromosome.isdigit()
                assert band.isdigit()
                #print (chromosome)
                #print (band)
                #print (subband)
                # check band and subband to see whether they're valid
            if chromosome.isdigit() and band.isdigit():
                chromosome = int(chromosome)
                #print(chromosome)
                if (1< chromosome <= 22) or (chromosome == "X" or "Y"):
                    return "True"
            else:
                return "False"
        else:
            return "False"
    else: 
        return "False"

is_valid_humanlocus("6p21.3")
is_valid_humanlocus("11q1.4")
is_valid_humanlocus("22p11.2")

is_valid_humanlocus("chr1:1000")
is_valid_humanlocus("nonsense")
is_valid_humanlocus("2a11p")

mystring = "21p11.2"
is_valid_humanlocus(mystring)

is_valid_humanlocus("6213.3")


# Part 2, a: Parse out the SNP id, chromosome, position and SNPs for each row.
# e.g. 

# snp = 'rs3094315chr1-742429(-,-)'
#line = snp.rstrip(')')
#line = snp.split('chr')
#snpID = line[0]
#line1 = line[1].split('(')
#line2 = line1[0].split("-")
#chromosome = line2[0]
#position = line2[1]
#snps = line1[1].strip(")")
#snps = snps.split(',')
#snp1 = snps[0]
#snp2 = snps[1]


# Open the sequence file in read mode.
#infile1 = open("GregMendel_SNPs.txt", 'r')
#infile2 = open("LillyMendel_SNPs.txt", 'r')

# create empty lists
#gchromosome = []
#gsnpID = []
#gposition = []
#gsnp1 =[]
#gsnp2 = []

#lchromosome = []
#lsnpID = []
#lposition = []
#lsnp1 = []
#lsnp2 = []

# loop to parse out the SNP information for Greg
#for line in infile1:
#    line = line.split('chr')
#    snpID = line[0]
#    line1 = line[1].split('(')
#   line2 = line1[0].split("-")
#    chromosome = line2[0]
#    position = line2[1]
#    snps = line1[1].strip(")")
#    snps = snps.split(',')
#    snp1 = snps[0]
#    snp2 = snps[1]
#    gchromosome.append(chromosome)
#    gsnpID.append(snpID)
#    gposition.append(position)
#    gsnp1.append(snp1)
#    gsnp2.append(snp2)
#print("Done")

# loop to parse out the SNP information for Lilly
#for line in infile2:
#    line = line.split('chr')
#    snpID = line[0]
#    line1 = line[1].split('(')
#    line2 = line1[0].split("-")
#    chromosome = line2[0]
#    position = line2[1]
#    snps = line1[1].strip(")")
#    snps = snps.split(',')
#    snp1 = snps[0]
#    snp2 = snps[1]
#    lchromosome.append(chromosome)
#    lsnpID.append(snpID)
#    lposition.append(position)
#   lsnp1.append(snp1)
#    lsnp2.append(snp2)
# print("Done")

# Part 2, b: define a function called read_SNP_file, which you then call from your main 
# script to process both Greg and Lilly’s data
#dict = {"ID":["rs3094315","rs111111",...],
#        "chr":[1,2,...],
#        "position":[1000,1009,...],
#        "SNP":["AG","AA",...]}

def read_SNP_file(filename):
    ''' Reads a SNP file and returns lists and a dictionary with parallel lists for each individual '''
    infile = open(filename,'r')
    ID = []
    chr = []
    position = []
    snp1 = []
    snp2 = []
    for line in infile:
        line = line.rstrip(')')
        line = line.split('chr')
        snpID = line[0]
        line1 = line[1].split('(')
        line2 = line1[0].split("-")
        chromosome = line2[0]
        assert chromosome.isdigit()
        pos = line2[1]
        snps = line1[1].rstrip(")")
        snps = snps.split(',')
        snps1 = snps[0]
        snps2 = snps[1]
        chr.append(chromosome)
        ID.append(snpID)
        position.append(pos)
        snp1.append(snps1)
        snp2.append(snps2)
    infile.close()
    dict = { "ID":[ID],
                 "chr": [chr],
                 "position": [position],
                 "SNP1":[snp1],
                 "SNP2":[snp2]}
    return [ID, chr, position, snp1, snp2, dict]
    
Lilly = read_SNP_file("LillyMendel_SNPs.txt")    
Greg = read_SNP_file("GregMendel_SNPs.txt")

lID = Lilly[0]
lchr = Lilly[1]
lposition = Lilly[2]
lsnp1 = Lilly[3]
lsnp2 = Lilly[4]
ldict = Lilly[5]

gID = Greg[0]
gchr = Greg[1]
gposition = Greg[2]
gsnp1 = Greg[3]
gsnp2 = Greg[4]
gdict = Greg[5]

#results = Lilly[5]
    

# Part 2, c:On Chromosome 10, find the largest region of shared SNPs between Lilly and Greg.

shared = []

for i, val in enumerate(lchr):
    if val == "10" and lsnp1[i] + lsnp2[i] == gsnp1[i] + gsnp2[i]:
        shared.append(lposition[i])
    else:
        continue
print ("Done")

shared.sort()
shared[0]
shared[22517]

# Part 2, d: Load the SNP definitions into a data structure so that you can look up a 
# description given a SNP id and the bases.

infile3 = open("SNP_definitions.txt", 'r')
line1 = infile3.readline()

SNP_def ={}

for line in infile3:
    line = line.strip('')
    line = line.split('\t')
    SNPid = line[0]
    genotype = line[1]
    des = line[2]
    SNP_def[SNPid] = genotype, des
print("Done")

# Part 2, e: identify SNPs from the region between 22070000 and 22106000 on chromosome 9 
# suggest about Greg’s chance of heart disease? What about Lilly’s chance of heart disease?

results = []

for i, val in enumerate(lchr):
    position = lposition[i]
    if val == "9" and 22070000 < int(position) > 22106000:
        results.append(lID[i])
    else:
        continue
print ("Done")

Heartdisease = {}
   
for k,v in SNP_def.items():
    if k == results and "heart disease" in v:
        key = k
        value = v
        Heartdisease[key] = value
    else:
        continue
print (Heartdisease)

#  They both do not have an increase risk of heart attack.

# Part 2, f: Find a SNP locus that interests you at SNPedia.com. Describe what is known 
# about the locus. Also, check what the SNP status is in both Lilly and Greg. 
# What does the SNP suggest about their health?

# SNP ID: rs6544713
# what is known: new locus identified for the increased risk of CAD
# minor allele of SNP (T) is associated with increased LDL
# on Chromosome	2
# Position	43846742
# Gene	ABCG8


for i, val in enumerate(lID):
    if val == "rs6544713":
        print (lsnp1[i] + lsnp2[i])
    else: 
        continue
print ("Done")

for i, val in enumerate(gID):
    if val == "rs6544713":
        print (gsnp1[i] + gsnp2[i])
    else: 
        continue
print ("Done")

# Since both Lilly and Greg have the SNP (T), they do have an increased risk of CAD.









