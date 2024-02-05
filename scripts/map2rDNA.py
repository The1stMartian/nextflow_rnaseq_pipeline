# Script is designed to map fastq files to the rDNA genome using Bowtie2
# This can be done with downsampled fastq files to rapidly get a sense of the pct of reads that are rRNA

import os
############################################################################################################
def fastq2sample(fq_folder_name, logfilepath, gnm, ff):
	'''	  
	  - All fastq file pairs should be formatted like "samplename_R1.fq" and "samplename_R2.fq"
	  - Collects all file names in fastq folder. 
	  - Identifies name of the sample from each fastq file pairs
	  - Returns list of sample names. These will be used for mapping.
	'''

	# collect sample names from the fastq file pairs in a folder
	list_of_names = []
	list_check = []
	fastaFileExtensions = []
	try:
		fnames = os.listdir(fq_folder_name)
		for i in fnames:
			i1 = i.split(".")[0].replace("_R1","")		# takes everything before '.' and removes _R1 or _R2
			i2 = i1.replace("_R2","")

			# get file extensions
			extension = i.split(".")[-1]
			if extension == "gz":
				extension = "." + i.split(".")[-2] + "." + i.split(".")[-1]

			if i2 not in list_check:
				list_check.append(i2) # prevents duplication
				list_of_names.append((i2, extension))
			fastaFileExtensions.append(extension)
		
		extensions = list(set(fastaFileExtensions))
		print("Processing the following samples:")
		for z in list_of_names:
			print(z[0])
		print("")

		# Record sample names in log
		with open(logfilepath, "w") as logf:
			logf.write("Genome Used: " + gnm + "\n")
			logf.write("Fasta File: " + ff + "\n")
			logf.write("Samples:\n")
			for z in list_of_names:
				logf.write(z[0] + "\n")
			logf.write("\n\n")
			logf.write("Command Log:\n")

		return list_of_names, extensions[0]

	except FileNotFoundError:
		print("Error: " + fq_folder_name + " not found. Please check the folder names and try again.")
		print()
		exit()

##########################################################################################################
def fastqcheck(sampleList, fqlocation):
	"""
	  - Checks that the fastq filenames are formatted as expected (X_R1.fq, and X_R2.fq) 
	  - Confirms that both forward and reverse files are there. 
	  - If files are not formatted correctly, script exits.
	  - If successful returns the file names.
	"""
	for s in sampleList:
		samplename = s[0]
		fqsuffix = s[1]
		f1 = samplename + "_R1." + fqsuffix
		f2 = samplename + "_R2." + fqsuffix
		filenames1 = os.listdir(fqlocation)

		if f1 in filenames1:
			if f2 in filenames1:
				continue
			else:
				print("Fastq files are not formatted properly for sample " + samplename + ".")
				print("Exiting.")
				exit()
	print("All fastq files detected in the input folder: " + fqlocation)
	print()
	return None

##########################################################################################################
def bowtie2(inputFilePathList, pwd, genome, fastqFldr, outFolder, logFilePath):
	"""
	- maps each fastq to the genome with Bowtie2
	- records mapping rate to file
	"""
	for sample in inputFilePathList:
		print("Mapping Sample: ", sample[0])
		f1 = pf + fastqFldr + sample[0] + "_R1." + sample[1]
		f2 = pf + fastqFldr + sample[0] + "_R2." + sample[1]
		c = "bowtie2 -x " + genome + " -1 " + f1 + " -2 " + f2 + " 2>> " + logFilePath + " |samtools view -bS - > " + outFolder + "/" + sample[0] + ".bam"
		print(c)
		with open(logFilePath, "a") as log:
			log.write("Sample: " + sample[0] + "\n")
			log.write("Files: " + sample[0] + "_R1." + sample[1] + ", " + sample[0] + "_R2." + sample[1] + "\n")
			log.write(c + "\n\n")
		os.system(c)
		with open(logFilePath, "a") as log:
			log.write("\n")
	return None

##########################################################################################################
def condenseLogFile(logFilePath, condensedPath):
	"""
	- Condenses log file into mini log file of sample:mapping rate
	"""
	outSamples = []
	outData = []
	with open(logFilePath, "r") as inFile:
		for line in inFile:
			if "Sample" in line:
				if "Samples" not in line:
					line = line.strip()
					outSamples.append(line.replace("Sample: ", "")) 
			elif "overall alignment rate" in line:
				line = line.replace("overall alignment rate\n", "").strip()
				outData.append(line)
	with open(condensedPath, "w") as outCondensed:
		outCondensed.write("Sample: ,")
		for l in outSamples:
			outCondensed.write(l + ",")
		outCondensed.write("\n")
		outCondensed.write("Overall Alignment Rate (rDNA): ,")
		for l in outData:
			outCondensed.write(l + ",")
		
	return None

#####################################
# Execute
#####################################
g = "YOUR_RDNA_GENOME_NAME"			# Name used in STAR
ff = "fastqSmall/" 					# Folder containing fasta files
fasta_file = "YOURFASTAFILE"		# Your rDNA fasta file name
genomeFileFldr = "genomefiles"		# Location of your fasta and STAR genome files
pf = os.getcwd() + "/"				# Current working directory
bamFolder = pf + "mapped2rDNA"		# BAM output folder

genomePath = pf + genomeFileFldr + "/" + g 

# Create outfolder
if not os.path.exists(bamFolder):
	os.mkdir(bamFolder)

# Create a command log file
logFilePath = pf + "logFile_rDNA.txt"
condensedLogFilePath = pf + "LogFile_rDNA_mini.csv"

print()
print("Expecting all paired-end files...")
print("Working directory is: ", pf)
print()
print("Mapping fastq files in " + ff + "...")
print()
listOfFqFiles, ext = fastq2sample(ff, logFilePath, g, fasta_file)
fastqcheck(listOfFqFiles, ff)
bowtie2(listOfFqFiles, pf, genomePath, ff, bamFolder, logFilePath)
condenseLogFile(logFilePath, condensedLogFilePath)