from Bio import SeqIO
import re
import sys

#Create empty lists to store amino acids in the first and second positions
x1 = []
x2 = []
x3 = []
x4 = []
masterlist = []

for e7 in SeqIO.parse("E7_aa_alpha_aligned.fasta", "fasta"):
        name = list(e7.id.split("_"))[0]
        sequence = e7.seq
        sequence2 = str(sequence)[82:95]
        sequence3 = ''.join(i for i in sequence2 if not i == '-')
        sequence4 = str(sequence)[117:123]
        match1 = re.search("C..C", str(sequence3))
        match2 = re.search("C..C", str(sequence4))
        if match1:
                if match2:
                        first = match1.group()[1]
                        second = match1.group()[2]
                        third = match2.group()[1]
                        fourth = match2.group()[2]
                        risk = list(e7.id.split("_"))[1]
                        x1.append(first)
                        x2.append(second)
                        info = [name, risk, first, second, third, fourth]
                        masterlist.append(info)
                        print("Name: " + name)
                        print("Risk: " + risk)
                        print("Motif: " + match1.group()[0:4] + "..." + match2.group()[0:4])
                        print("X1: " + first)
                        print("X2: " + second)
                        print("X3: " + third)
                        print("X4: " + fourth + "\n")
                else:
                        print("No match for second cxxc in " + name + "\n")
        else:
                print("No match for first cxxc in " + name + "\n")


highlist = [entry for entry in masterlist if entry[1] == "H"]
x1high = []
x2high = []
x3high = []
x4high = []
for entry in highlist:
        x1high.append(entry[2])
        x2high.append(entry[3])
        x3high.append(entry[4])
        x4high.append(entry[5])
bothlist = [entry for entry in masterlist if entry[1] == "B"]
x1both = []
x2both = []
x3both = []
x4both = []
for entry in bothlist:
        x1both.append(entry[2])
        x2both.append(entry[3])
        x3both.append(entry[4])
        x4both.append(entry[5])
lowlist = [entry for entry in masterlist if entry[1] == "L"]
x1low = []
x2low = []
x3low = []
x4low = []
for entry in lowlist:
        x1low.append(entry[2])
        x2low.append(entry[3])
        x3low.append(entry[4])
        x4low.append(entry[5])

with open("cxxc_alpha_output.txt", "w") as t:
        sys.stdout = t
        print("onco")
        for entry in highlist:
                motif = ("C",entry[2],entry[3],"C-C",entry[4],entry[5],"C")
                print(''.join(str(x) for x in motif))
        print("possibly")
        for entry in bothlist:
                motif = ("C",entry[2],entry[3],"C-C",entry[4],entry[5],"C")
                print(''.join(str(x) for x in motif))
        print("low")
        for entry in lowlist:
                motif = ("C",entry[2],entry[3],"C-C",entry[4],entry[5],"C")
                print(''.join(str(x) for x in motif))

input()
