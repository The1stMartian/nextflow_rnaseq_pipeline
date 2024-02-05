# Script can be used to format a fasta sequence with numbers and spaces
# Script will remove all numbers and spaces

outlines = []
removelist = ["1", "2", "3", "4", "5", "6", "7", "8", "9","0", " "]
with open("INPUTFILEPATH", "r") as f:
	for i in f:
		newline = ""
		for l in i:
			if l not in removelist:
				newline = newline + l
		outlines.append(newline)

with open("OUTPUTFILEPATH", "w") as out:
	for n in outlines:
		out.write(n)
exit()
