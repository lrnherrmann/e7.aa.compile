from Bio import SeqIO
import re
import sys

#Create empty lists to store amino acids in the first second and third positions
x1 = []
x2 = []
x3 = []
masterlist = []

for e7 in SeqIO.parse("E7_DNA_alpha.fasta", "fasta"):
        name = list(e7.id.split("|"))[3]
        sequence = e7.seq
        translation = sequence.translate()
        match = re.search(".L.C.E", str(translation))
        if match:
                first = match.group()[0]
                second = match.group()[2]
                third = match.group()[4]
                risk = list(e7.id.split("|"))[-1]
                x1.append(first)
                x2.append(second)
                x3.append(third)
                info = [name, risk, first, second, third]
                masterlist.append(info)
                print("Name: " + name)
                print("Risk: " + risk)
                print("Motif: " + match.group()[0:6])
                print("X1: " + first)
                print("X2: " + second)
                print("X3: " + third)
        else:
                print("No match for " + name + "\n")


highlist = [entry for entry in masterlist if entry[1] == "H"]
x1high = []
x2high = []
x3high = []
for entry in highlist:
        x1high.append(entry[2])
        x2high.append(entry[3])
        x3high.append(entry[4])
bothlist = [entry for entry in masterlist if entry[1] == "B"]
x1both = []
x2both = []
x3both = []
for entry in bothlist:
        x1both.append(entry[2])
        x2both.append(entry[3])
        x3both.append(entry[4])
lowlist = [entry for entry in masterlist if entry[1] == "L"]
x1low = []
x2low = []
x3low = []
for entry in lowlist:
        x1low.append(entry[2])
        x2low.append(entry[3])
        x3low.append(entry[4])

with open("XLXCXE_alpha_output.txt", "w") as t:
        sys.stdout = t
        for entry in masterlist:
                print(entry[0],"(",entry[1],"): ",entry[2]," ",entry[3]," ",entry[4])
                
        print("\nX1 amino acid frequencies:")
        print("\nOverall:")
        for aa in set(x1):
                percent = 100*(x1.count(aa) / len(x1))
                print(aa, ": ", x1.count(aa), "(", round(percent, 2), "%)")
        print("\nHigh risk:")
        for aa in set(x1high):
                percent = 100*(x1high.count(aa) / len(x1high))
                print(aa, ": ", x1high.count(aa), "(", round(percent, 2), "%)")
        print("\nLimited evidence:")
        for aa in set(x1both):
                percent = 100*(x1both.count(aa) / len(x1both))
                print(aa, ": ", x1both.count(aa), "(", round(percent, 2), "%)")
        print("\nLow risk:")
        for aa in set(x1low):
                percent = 100*(x1low.count(aa) / len(x1low))
                print(aa, ": ", x1low.count(aa), "(", round(percent, 2), "%)")

        print("\n\nX2 amino acid frequencies:")
        print("\nOverall:")
        for aa in set(x2):
                percent = 100*(x2.count(aa) / len(x2))
                print(aa, ": ", x2.count(aa), "(", round(percent, 2), "%)")
        print("\nHigh risk:")
        for aa in set(x2high):
                percent = 100*(x2high.count(aa) / len(x2high))
                print(aa, ": ", x2high.count(aa), "(", round(percent, 2), "%)")
        print("\nLimited evidence:")
        for aa in set(x2both):
                percent = 100*(x2both.count(aa) / len(x2both))
                print(aa, ": ", x2both.count(aa), "(", round(percent, 2), "%)")
        print("\nLow risk:")
        for aa in set(x2low):
                percent = 100*(x2low.count(aa) / len(x2low))
                print(aa, ": ", x2low.count(aa), "(", round(percent, 2), "%)")
                
        print("\n\nX3 amino acid frequencies:")
        print("\nOverall:")
        for aa in set(x3):
                percent = 100*(x3.count(aa) / len(x3))
                print(aa, ": ", x3.count(aa), "(", round(percent, 2), "%)")
        print("\nHigh risk:")
        for aa in set(x3high):
                percent = 100*(x3high.count(aa) / len(x3high))
                print(aa, ": ", x3high.count(aa), "(", round(percent, 2), "%)")
        print("\nLimited evidence:")
        for aa in set(x3both):
                percent = 100*(x3both.count(aa) / len(x3both))
                print(aa, ": ", x3both.count(aa), "(", round(percent, 2), "%)")
        print("\nLow risk:")
        for aa in set(x3low):
                percent = 100*(x3low.count(aa) / len(x3low))
                print(aa, ": ", x3low.count(aa), "(", round(percent, 2), "%)")


input()
