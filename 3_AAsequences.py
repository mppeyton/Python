#!/usr/bin/python
"""
@author: minapeyton

"""

# Script: charged_sequence.py
# Problem Description: 
#    Read in amino acid sequences from an input file and
#    calculate statistics on the 'disordered' region.


# Open the sequence file in read mode.
infile = open("Hsp90_conserved.txt", 'r')

# Open an empty file called 'aaoutput.txt' in write mode
outfile = open("aaoutput.txt", 'wt')

##### First, develop the code to process just 1 sequence before trying to
##### process all lines in the file. Use the .readline() method to read the first
##### line of the file and perform all operations and calculations only
##### on that line. Once your code for the first line works, paste it into the 
##### Section 2 below to make it run on all lines.


##### Code for processing the first line only
# Read in the first line of the file
line1 = infile.readline()

# Remove any trailing newline characters from the sequence
line1 = line1.rstrip()

# Split line on tab character and save into two variables
fragments = line1.split(sep='\t') # list of 2 fragments

# fragment 1 is the name
name = ""
name = fragments[0]
name = name.rstrip()

# fragment 2 is the sequence
seq = ""
seq = fragments[1]
seq = seq.rstrip()

# Extract the 'disordered' region of the Hsp90 sequence

# first sequence fragment split at EYLEE at the start of disordered region
seqfrag1 = seq.split('EYLEE')
seqfrag1

# sequence with disordered region
seqfrag2 = seqfrag1[1]

# split disordered region out of seqfrag2 at LNKTKP
seqfrag2 = seqfrag2.split('LNKTKP')
seqfrag2

# extract disordered region
disordered = seqfrag2[0]
disordered


# Find the percentage of negatively charged amino acids in the disordered region
# there are two negatively charged amino acids: aspartic acide (Asp, D) and 
# glutamic acid (Glu, E)

# count negatively charged amino acid Asp (D)
Asp = disordered.count('D') #11

# count negatively charged amino acid Glu (E)
Glu = disordered.count('E') #21

# percent of Asp and Glu in the disordered region
# length of the disordered region
disorderedlen = len(disordered) #86

Asppercent = ((Asp/disorderedlen)*100) #12.79%
Glupercent = ((Glu/disorderedlen)*100) #24.42%

# Find the percentage of positively charged amino acids in the disordered region
# there are three positively charged amino acids: lysine (Lys, K), arginine (Arg, R)
# and histidine (His, H)

# count positively charged amino acid Lys (K)
Lys = disordered.count('K') #19

# count positively charged amino acid Arg (R)
Arg = disordered.count('R') #3

# count positively charged amino acid His (H)
His = disordered.count('H') #1

# percent of Lys, Arg, and His in the disordered region

Lyspercent = ((Lys/disorderedlen)*100) #22.09%
Argpercent = ((Arg/disorderedlen)*100) #3.49%
Hispercent = ((His/disorderedlen)*100) #1.16%


# Repeat the calculations above for the complete sequence

# length of the complete sequence
seqlen = len(seq) #724

# count of negatively charged amino acid in complete sequence
seqAsp = seq.count('D') #51
seqGlu = seq.count('E') #96

# count of positively charged amino acid in complete sequence
seqLys = seq.count('K') #75
seqArg = seq.count('R') #32
seqHis = seq.count('H') #13

# percent of negatively charged amino acid in complete sequence
seqAsppercent = ((seqAsp/seqlen)*100) #7.04%
seqGlupercent = ((seqGlu/seqlen)*100) #13.26%

# percent of positvely charged amino acid in complete sequence

seqLyspercent = ((seqLys/seqlen)*100) #10.36%
seqArgpercent = ((seqArg/seqlen)*100) #4.42%
seqHispercent = ((seqHis/seqlen)*100) #1.80%

# print output for this line to the file

# fraction of disordered negatively charged amino acids
frac_disNeg = Asppercent + Glupercent #37.21%

# fraction of disordered positively charged amino acids
frac_disPos = Lyspercent + Argpercent + Hispercent #26.74%

# fraction of sequence negatively charged amino acids
frac_seqNeg = seqAsppercent + seqGlupercent #20.30%

# fraction of sequence positively charged amino acids
frac_seqPos = seqLyspercent + seqArgpercent + seqHispercent #16.57

# convert all fraction objects into a string variable instead of a float 
# in order to print

frac_disPos = str(frac_disPos)
frac_disNeg = str(frac_disNeg)

disorderedlen = str(disorderedlen)

frac_seqPos = str(frac_seqPos)
frac_seqNeg = str(frac_seqNeg)


# print out:
# The name of the amino acid sequence taken from file.
# Percentage of positively charged amino acids in the disordered region. 
# Percentage of negatively charged amino acids in the disordered region. 
# Length of disordered region.
# Percentage of positively charged amino acids in the whole sequence.
# Percentage of negatively charged amino acids in the whole sequence.

print("name", "frac_disPos", 'frac_disNeg', 'disorderedlen', 'frac_seqPos', 'frac_seqNeg',sep='\t',end='\n', file = outfile)
print(name, frac_disPos, frac_disNeg, disorderedlen, frac_seqPos, frac_seqNeg,sep='\t',end='\n', file = outfile)

# close infile and outfile
infile.close()
outfile.close()
# must close out outfile in order for the file content to save


##### Section 2:  Code for processing all lines in the file
##### NOTE: only try this after you've gotten your code for the
##### first line to work above.
##### To get the for loop to work, replace the XXX and YYY with the appropriate variable names,
##### then copy your code from above. You'll need to uncomment the line with "for".

# Open the sequence file in read mode.
infile = open("Hsp90_conserved.txt", 'r')

# Open an empty file called 'aaoutput.txt' in write mode
outfile = open("aaoutput.txt", 'wt')
print("name", "frac_disPos", 'frac_disNeg', 'disorderedlen', 'frac_seqPos', 'frac_seqNeg',sep='\t',end='\n', file = outfile)

for line in infile:
    line = line.rstrip()
    fragments = line.split(sep='\t')
    name = ''
    name = fragments[0]
    name = name.rstrip()
    seq = ''
    seq = fragments[1]
    seq = seq.rstrip()
    if seq.find("EYLEE" and "LNKTKP"):
        seqfrag1 = seq.split('EYLEE')
        seqfrag2 = seqfrag1[1]
        seqfrag2 = seqfrag2.split('LNKTKP')
        disordered = seqfrag2[0]
    else:
        continue
    disorderedlen = len(disordered)
    # count negatively charged amino acid Asp (D)
    Asp = disordered.count('D')
    # count negatively charged amino acid Glu (E)
    Glu = disordered.count('E')
    # percent of Asp and Glu in the disordered region
    Asppercent = ((Asp/disorderedlen)*100) 
    Glupercent = ((Glu/disorderedlen)*100) 
    # count positively charged amino acid Lys (K)
    Lys = disordered.count('K') 
    # count positively charged amino acid Arg (R)
    Arg = disordered.count('R') 
    # count positively charged amino acid His (H)
    His = disordered.count('H') 
    # percent of Lys, Arg, and His in the disordered region
    Lyspercent = ((Lys/disorderedlen)*100) 
    Argpercent = ((Arg/disorderedlen)*100) 
    Hispercent = ((His/disorderedlen)*100) 
    # length of the complete sequence
    seqlen = len(seq) #724
    # count of negatively charged amino acid in complete sequence
    seqAsp = seq.count('D') #51
    seqGlu = seq.count('E') #96
    # count of positively charged amino acid in complete sequence
    seqLys = seq.count('K') #75
    seqArg = seq.count('R') #32
    seqHis = seq.count('H') #13
    # percent of negatively charged amino acid in complete sequence
    seqAsppercent = ((seqAsp/seqlen)*100) #7.04%
    seqGlupercent = ((seqGlu/seqlen)*100) #13.26%
    # percent of positvely charged amino acid in complete sequence
    seqLyspercent = ((seqLys/seqlen)*100) #10.36%
    seqArgpercent = ((seqArg/seqlen)*100) #4.42%
    seqHispercent = ((seqHis/seqlen)*100) #1.80%
    # fraction of disordered negatively charged amino acids
    frac_disNeg = Asppercent + Glupercent #37.21%
    # fraction of disordered positively charged amino acids
    frac_disPos = Lyspercent + Argpercent + Hispercent #26.74%
    # fraction of sequence negatively charged amino acids
    frac_seqNeg = seqAsppercent + seqGlupercent #20.30%
    # fraction of sequence positively charged amino acids
    frac_seqPos = seqLyspercent + seqArgpercent + seqHispercent #16.57
    # convert all fraction objects into a string variable instead of a float 
    # in order to print
    frac_disPos = str(frac_disPos)
    frac_disNeg = str(frac_disNeg)
    disorderedlen = str(disorderedlen)
    frac_seqPos = str(frac_seqPos)
    frac_seqNeg = str(frac_seqNeg)
    print(name[1:], frac_disPos[0:4], frac_disNeg[0:4], disorderedlen, frac_seqPos[0:4], frac_seqNeg[0:4],sep='\t',end='\n', file = outfile)

print("Done")

# Close the input file -- (NOTE this is not inside the for loop, so it should not be indented)
infile.close()
# Close the output file
outfile.close()
