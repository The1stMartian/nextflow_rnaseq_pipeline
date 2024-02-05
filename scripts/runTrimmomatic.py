# An easy-to-use wrapper for Trimmomatic

import os

print()
print("Trim fastq files with Trimmomatic:")
print("  - This script assumes all file pairs are named like 'sample1_R1' and 'sample1_R2'.")
print("  - Fastq files can have extensions  '.fq', or '.fastq' and can be gzipped or not (ending .fastq.gz or .fq.gz)")
print("  - The script is also expecting all fastq files to have the same file extensions. It's not set up for a mix.")
print()

print("Where are the fastq files that need to be trimmed?")
inputPath = input("Please provide the relative path (e.g. 'fastq'): ")
outputPath = input("What is the output folder's relative path? ")

print("Ok, the original fastq files are in ", inputPath)
print("and the output files will be saved to: ", outputPath)

def getFastqFiles(inputFolder, F1ext, F2ext):
	"""
	- From input folder, get all fastq files
	"""
	fastqExtList = [".fastq", ".fq", ".fastq.gz", "fq.gz"]
	try:
		fnames = os.listdir(inputFolder)
		sampleList = []
		endingList = []
		for f in fnames:
			splt = f.split(".")
			end1 = "." + splt[-1]
			end2 = "." + splt[-2] + "." + splt[-1]
			if end1 in fastqExtList:
	
				# Remove _R1 and _R2 from sample names to get generic sample name
				sampleName = sampleName.replace(F1ext, "")
				sampleName = sampleName.replace(F2ext, "")
				if sampleName not in sampleList:
					sampleList.append(sampleName)

				if end1 not in endingList:
					endingList.append(end1)

			elif end2 in fastqExtList:
				# Remove _R1 and _R2 from sample names to get generic sample name
				sampleName = f.split(".")[0]
				
				# Remove _R1 and _R2 from sample names to get generic sample name
				sampleName = sampleName.replace(F1ext, "")
				sampleName = sampleName.replace(F2ext, "")
				if sampleName not in sampleList:
					sampleList.append(sampleName)
				
				if end2 not in endingList:
					endingList.append(end2)

			else:
				continue
		if len(endingList) == 1:
			print()
			print("All detected files end with: ", endingList[0])
			print("Samples:")
			for f in sampleList:
				print(f)
			print()
		elif len(endingList) != 1:
			print("Multiple file endings detected:")
			for e in endingList:
				print(e)
			print("Please ensure that files are named consistently. Only one file type is expected.")
			exit()

		return sampleList, endingList[0]
	except FileNotFoundError:
		print("Input folder not found. Exiting...")
		print()
		exit()

def runTrimmomatic(fileList, wd, inputFolder, outputFolder, fileExt, ending1, ending2, logFileName, scriptsFldr, adapterFile):
	"""
	- INput is a list of fastq files, the relative input folder path, and the present working directory
	  which should be the parent directory of the fastq folder
	- runs Trimmomatic on a list of fastq files, saves to output folder input by user
	- removes truseq universal adaptors in fastq files
	"""
	# Make the output (trimmed fastq) folder
	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)

	# Trimmomatic Paths
	trimmomaticPath = wd + scriptsFldr + "/" + "trimmomatic-0.39.jar"
	adapterFilePath = wd + scriptsFldr + "/" + "adapters/" + adapterFile

	print()
	print("Expecting to find Trimmomatic at location:")
	print("   " + trimmomaticPath)
	print()
	print("Expecting to find adapter sequences at:")
	print("   " + adapterFilePath)

	# Trim files and record commands
	for f in fileList:
		
		# Input files
		forwardReadFile = inputFolder + f + ending1 + fileExt
		reverseReadFile = inputFolder + f + ending2 + fileExt
		
		# Output files
		output1paired = outputFolder + f + "_paired" + ending1 + fileExt
		output2paired = outputFolder + f + "_paired" + ending2 + fileExt
		output1unpaired = outputFolder + f + "_unpaired" + ending1 + fileExt
		output2unpaired = outputFolder + f + "_unpaired" + ending2 + fileExt
		
		print("Opening input files:")
		print("  " + forwardReadFile)
		print("  " + reverseReadFile)
		print()

		print("Outputting files:")
		print("  " + output1paired)
		print("  " + output2paired)
		
		c = "java -jar " + trimmomaticPath + " PE "
		c = c + forwardReadFile + " " + reverseReadFile + " "
		c = c + output1paired	+ " " + output1unpaired + " "
		c = c + output2paired	+ " " + output2unpaired + " "
		c = c + "ILLUMINACLIP:" + adapterFilePath + ":2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 2>> " + logFileName
		print()
		print("Issuing the following command to the command line:")
		print("  " + c)
		os.system(c)
	print()
	print("Trimming complete.")
	print()

##################################################
# Execute
##################################################
pwd = os.getcwd() + "/"

# Clean input folder name:
if inputPath.startswith("./"):
	inputPath = inputPath[2:]

# Input info
fastqInputDirectory = pwd + inputPath + "/"
trimmedDirectory = pwd + outputPath + "/"
file1ext = "_R1"
file2ext = "_R2"
logfilename = "./Log_Trimmomatic.txt"
scripts = "scripts"
adapterfile = "TruSeq3-PE-2.fa"

# Run...
fqFileList, fqFileExtension = getFastqFiles(fastqInputDirectory, file1ext, file2ext)
runTrimmomatic(fqFileList, pwd, fastqInputDirectory, trimmedDirectory, fqFileExtension, file1ext, file2ext, logfilename, scripts, adapterfile)