# Script copies the first N number of reads to a new file
# The smaller files can be used to quickly map to rDNA (or other) to test
# for rRNA (or tRNA, etc.) contamination

import os
import csv

def getpaths(inputDir):
	"""
	- From input directory name, collect all .fastq file paths
	- Shrink fastq files to output read number
	- Print smaller fastq files into output folder
	"""
	filepaths = []
	fnames = os.listdir(inputDir)
	for f in fnames:
		splt = f.split(".")
		
		gzCheck = (splt[-2] + "." + splt[-1])
		if gzCheck == "fastq.gz":
			fileExt = ".fastq.gz"
		elif gzCheck == "fq.gz":
			fileExt = ".fq.gz"
		
		else:
			if splt[-1] == "fastq":
				fileExt = ".fastq"
			elif splt[-1] == "fq":
				fileExt = ".fq"
		
		inputFilePath = inputDir + "/" + f
		filepaths.append((f, inputFilePath, fileExt))
	return filepaths


def shrinkFastas(pathList, outReadNumber, outputDir):
	"""
	- Shrinks files down to specified number of reads
	- If inputs are .gz, so are outputs
	"""
	uncompressedList = [".fastq", ".fq"]
	compressedList = [".fastq.gz", ".fq.gz"]
	totalLines = outReadNumber*4
	for f in pathList:
		print("Working on: ", f[0])
		
		# For non-compressed fastas
		if f[2] in uncompressedList:
			outFilePath = outputDir + "/" + f[0]
			outlines = []
			with(open(outFilePath, "w")) as out:
				with open(f[1], "r") as inFile:
					reader = csv.reader(inFile)
					for idx, row in enumerate(reader):
						if idx <= totalLines:
							outlines.append(row)
						else:
							break
				for l in outlines:
					out.write(l[0] + "\n")

		# For compressed fastas
		elif f[2] in compressedList:
			outFile = f[0].split(".")[0] + ".fastq"
			outFilePath = outputDir + "/" + outFile
			c = "gunzip -c " + f[1] +  " | head -n " + str(totalLines) + " > " + outFilePath
			print("command: ", c)
			os.system(c)

			print()
			print("Ignore any gunzip 'Broken Pipe' error. It works just fine. It's a python/os thing.")


#######################################################
# Execute
#######################################################
fqFolder = "./fastq"
print()
print("This script will shrink all fastq files found in the ", fqFolder)
finalNumberLines = int(input("How many reads do you want in the outfiles?"))
pathList = getpaths(fqFolder)
shrinkFastas(pathList, finalNumberLines, "./fastqSmall")
print("Script complete.")
print()