# Requiremes previous installation of the following in Windows Subsystem for Linux: 
# 1) STAR aligner (apt install rna-star)
# 2) A STAR genome of choice
# 3) A fasta formatted genome file and accompanying features file (.gtf)
# 4) A folder called "fastq" in the same directory as this script. Files should be named like "sample1_R1.fastq.gz" and "sample1_R2.fastq.gz"
#	 The change_filenames.py script is provided to assist with name formatting.
# 5) A folder named "scripts" in the same directory. Should contain the featureCounts executable and FastQC app folder.

# Run using with Python3 on linux without any arguments.
# Change the genome within the script

import os
import subprocess
import shutil
############################################################################################################
def fastq2sample(fq_folder_name):
	'''	  
	  - All fastq file pairs should be formatted like "samplename_R1.fq" and "samplename_R2.fq"
	  - Collects all file names in fastq folder. 
	  - Identifies name of the sample from each fastq file pairs
	  - Returns list of sample names. These will be used for mapping.
	'''

	# collect sample names from the fastq file pairs in a folder
	list_of_names = []
	fastaFileExtensions = []
	try:
		fnames = os.listdir(fq_folder_name)
		for i in fnames:
			i1 = i.split(".")[0].replace("_R1","")		# takes everything before '.' and removes _R1 or _R2
			i2 = i1.replace("_R2","")
			if i2 not in list_of_names:
				list_of_names.append(i2)

			# get file extensions
			extension = i.split(".")[-1]
			if extension == "gz":
				extension = "." + i.split(".")[-2] + "." + i.split(".")[-1]
			fastaFileExtensions.append(extension)
		extensions = list(set(fastaFileExtensions))
		print("Fastq file extensions: ")
		for e in extensions:
			print(e)
		if len(extensions) != 1:
			print("ERROR: Multiple Fastq file extensions detected! Please use consistent fastq file formats and repeat.")
			print()
			#exit()
		elif len(extensions) == 1:
			print("Fastq file extensions are all: ", extensions[0])
			print("Contining...")

		print("Processing the following samples:")
		logFile.write("Samples:\n")
		for z in list_of_names:
			print(z)
		print("")

		# Record sample names in log
		for z in list_of_names:
			logFile.write(z + "\n")
		logFile.write("\n")

		return list_of_names, extensions[0]

	except FileNotFoundError:
		print("Error: " + fq_folder_name + " not found. Please check the folder names and try again.")
		print()
		exit()

##########################################################################################################
def fastqcheck(samplename, fqlocation, fqsuffix):
	"""
	  - Checks that the fastq filenames are formatted as expected (X_R1.fq, and X_R2.fq) 
	  - Confirms that both forward and reverse files are there. 
	  - If files are not formatted correctly, script exits.
	  - If successful returns the file names.
	"""

	f1 = samplename + "_R1" + fqsuffix
	f2 = samplename + "_R2" + fqsuffix
	filenames1 = os.listdir(fqlocation)

	if f1 in filenames1:
		if f2 in filenames1:
			return f1, f2
		else:
			print("Fastq files are not formatted properly for sample " + samplename + ".")
			print("Exiting.")
			exit()
	return None
	
##########################################################################################################
def filecheck(rc_fldr, mapf, saf_file, gfldr, fasta_file, scripts_fldr, fq_fldr, fqcYN, mapYN, fcYN):
	'''
	  - Makes sure all of the expected folders and files are there before mapping
	  	- rc_fldr = readcounts output folder
	  	- gflder = genomes file folder
	  	- .saf file is a tab delimited gene coordinates file. Users can create their own using
	  	  existing .saf files by copy/pasting the file and replacing the data colums in Excel.
	'''
	# Checks on existence of folders/files/etc
	if not os.path.exists(mapf):		# make mapped folder
		os.makedirs(mapf)
	if not os.path.exists(rc_fldr):		# make readcounts folder
		os.mkdir(rc_fldr)
	if not os.path.exists(scripts_fldr):
		print("Script folder is missing. Replace and restart.")
		exit()
	if not os.path.exists(gfldr + saf_file):
		print("The .saf/.gtf file is missing or improperly named. Add it to the genome files folder and re-try.")
		print("Looking for: ", gfldr + saf_file)
		exit()
	
	# Check for .fasta file
	# Possible fasta file names
	fasta_file_name = fasta_file + ".fasta"
	alt_fasta_file1 = fasta_file + ".fna"
	alt_fasta_file2 = fasta_file + ".fa" 

	if not os.path.exists(gfldr + fasta_file_name):
		if os.path.exists(gfldr + alt_fasta_file1):
			fasta_file_name = alt_fasta_file1
		else:
			if os.path.exists(gfldr + alt_fasta_file2):
				fasta_file_name = alt_fasta_file1
			else:
				print("The .fasta (or .fna) file is missing or improperly named. Add it to the genome files folder and re-try.")
				exit()
	
	# If FastQC or Mapping is "yes", look for fastq files:
	if fqcYN in yesList or mapYN in yesList:
		# Gets sample names from fastq file pairs in the fastq folder
		fq_files, fqsuffix = fastq2sample(fq_fldr)
		fq_file_names = []
		for name in fq_files:

			# Confirms expected formatting of fastq files from a given pair
			fq1, fq2 = fastqcheck(name, fq_fldr, fqsuffix)
			
			if fq1 != "error":
				if fq2 != "error":
					fq_file_names.append([fq1,fq2,name.strip("_")])		# Stores names as [name_1.fq, name_2.fq, name]
				else:
					print("Warning! Cannot find fastq file #2: ", fq2)
					print("Check that the sample names are accurate, and that")
					print("files are have '..._R1.fq' or '..._R1.fastq'")
					exit()
			else:
				print("Cannot find fastq file #1: " + fq1)
				print("Check that the sample names are accurate, and that")
				print("files are have '..._R1.fq' or '..._R1.fastq'")
				exit()

		print()
		print("File check completed successfully!")
		print("Continuing...")
		print()
		return fq_file_names, fasta_file_name

##########################################################################################################
def fastqc(fqList, fastqcFldr, pwd, fqFolder, logfilepath):
	"""
	- Runs fastqc on the input fastq files
	- FastQC can run in parallel up to the number of files specified by the -t flag (1 file per thread)
	- Flags --noextract prevents extraction of out files, 
	"""
	print("Running FastQC: ")
	
	# Initial command
	c = "fastqc --noextract -t 6 -o " + pwd + fastqcFldr 

	# Add all files paths, serially, to the command
	for f in fqList:
		c = c + " " + pwd + fqFolder + f[0] + " " + pwd + fqFolder + f[1] + " >> " + logfilepath
		print("c")
	print("CMD: ", c)
	logFile.write(c + "\n")
	os.system(c)
	logFile.write("\n")

	return None

##########################################################################################################3
def mapFQ(file_pair, fq_folder, star_gnm_fldr, pwd, starOutFldr, tmp, logfilepath):
	"""
	  - Runs fastqc on all fastq files
	  - Maps paired-end fastq files to the indicated genome (bottom) by issuing commands to the terminal.
	  - For each pair of paired-end fastq files in the fastq folder:
		- Maps both files to the user-specified genome (below) using STAR
	  	- Uses featureCounts > read counts (quantitative table of feature:mapped read number) per feature
	"""

	# Zipped fastq file names
	fqf1 = file_pair[0]		# forward fastq file
	fqf2 = file_pair[1]		# reverse fastq file
	name = file_pair[2]		# generic file name

	# Fastq file names uncompressed (without .gz)
	fqfa = fqf1.replace(".gz","")
	fqfb = fqf2.replace(".gz","")

	# Unzip before running b/c STAR is being finicky
	print("Unzipping fastq files:")
	c = "gunzip -c " + pwd + fq_folder + fqf1 + " > " + pwd + fq_folder + fqfa
	print(c)
	os.system(c)
	c = "gunzip -c " + pwd + fq_folder + fqf2 + " > " + pwd + fq_folder + fqfb
	print(c)
	os.system(c)

	print("Using STAR to map: " + name)

	# Ensure Tmp folder is not present or STAR throws an error
	if os.path.exists(tmp):
		shutil.rmtree(tmp)

	# STAR mapping parameters can be customized here
	c = "STAR --genomeDir " + \
		star_gnm_fldr + \
		" --runThreadN 6 " + \
		" --readFilesIn " + pwd + fq_folder + fqfa + " " + pwd + fq_folder + fqfb + \
		" --outSAMtype BAM SortedByCoordinate " + \
		" --outTmpDir " + tmp + \
		" --quantMode GeneCounts " + \
		" --outFileNamePrefix " + starOutFldr + name + "." + " >> " + logfilepath 

	# Run STAR
	print("CMD: ", c)
	logFile.write(c + "\n")
	os.system(c)
	logFile.write("\n")

	# Delete the uncompressed fastq files
	print("Deleting uncompressed fastq file: " + pwd + fq_folder + fqfa)
	logFile.write("Removing uncompressed fastq file via os.remove: " + pwd + fq_folder + fqfa)
	os.remove(pwd + fq_folder + fqfa)
	print("Deleting uncompressed fastq file: " + pwd + fq_folder + fqfb)
	logFile.write("Removing uncompressed fastq file via os.remove: " + pwd + fq_folder + fqfb)
	os.remove(pwd + fq_folder + fqfb)

	print("Mapping complete for sample: " + name)
	return None

##########################################################################################################
def featurecounts(mappedfolder, genomeFilesFldr, rcf, fFile):
	""" 
	- Runs featureCounts on all .bam files found in the mapped folder
	"""
	# Collect the .bam files output by STAR (everything in the mapped folder)
	onlyfiles = [f for f in os.listdir(mappedfolder) if os.path.isfile(os.path.join(mappedfolder, f))]
	bamfiles = [f for f in onlyfiles if ".bam" in f]

	# Feature counts is called here
	for f in bamfiles:
		outFileName = rcf + f.split(".")[0] + ".txt"
		c = "scripts/featureCounts -a " + genomeFilesFldr + fFile + " -F GTF -o " + outFileName + " " + mappedfolder + f 
		print("CMD: ", c)
		logFile.write(c + "\n")
		os.system(c)
	return None

##########################################################################################################
### User Customization
########################################################################################################## 
g = "dm6" 							# Genome name (the name of the STAR genome folder)
pf = os.getcwd() + "/"				# Parent folder, script location
ff = "fastq/" 						# Folder where trimmed/cleaned .fq files are
mf = "mapped/"						# Output mapping folder
rf = "read_counts/"					# Output folder for readcounts
so = "mapped/"						# Folder for mapped/aligned .bam/.sam files
gf = "genomefiles/"					# Genome .fasta and .gtf files
sf = "scripts/"						# Scripts/app folder
fqc = "fastqc/" 					# fastqc output folder name
tf = "temp/"						# Temp files folder
fasta_file = g 						# Name of the fasta formatted genome sequence file (can be .fna or .fasta)
featuresFile = g + ".gtf"			# Name of features for featureCounts
fastq_suffix = ".fastq.gz"			# Fastq file suffix
gf = pf + gf 						# Full genome files folder path
sg = pf + gf + g 					# Full STAR genome path
temp = pf + tf
##########################################################################################################
# Execute:
##########################################################################################################
print()
print(" ********                     **                 **                         **                  **        ")
print("/**/////                     /**    ****        ****                       /**  **   **        //  ")
print("/**        ******    ****** ****** **//**      **//**   *******   ******   /** //** **   ****** **  ******")
print("/*******  //////**  **//// ///**/ /** /**     **  //** //**///** //////**  /**  //***   **//// /** **//// ")
print("/**////    ******* //*****   /**  //*****    ********** /**  /**  *******  /**   /**   //***** /**//***** ")
print("/**       **////**  /////**  /**   ////**   /**//////** /**  /** **////**  /**   **     /////**/** /////**")
print("/**      //******** ******   //**     /***  /**     /** ***  /**//******** ***  **      ****** /** ****** ")
print("//        //////// //////     //      ///   //      // ///   //  //////// ///  //      //////  // //////  ")
print("")
print(" *******  **                  ** **   ")
print("/**////**//  ******          /**// ")
print("/**   /** **/**///**  *****  /** ** *******   *****  ")
print("/******* /**/**  /** **///** /**/**//**///** **///**  ")
print("/**////  /**/****** /******* /**/** /**  /**/*******  ")
print("/**      /**/**///  /**////  /**/** /**  /**/**////  ")
print("/**      /**/**     //****** ***/** ***  /**//****** ")
print("//       // //       ////// /// // ///   //  //////  ")
print("")

# Collect run instructions from the user
answer0 = input("Run FastQC on all fastq files? (Y/N): ")
answer1 = input("Should the fastq files be mapped to genome '" + g  + "'? (Y/N): ")
answer2 = input("Use featureCounts? (Y/N): ")
yesList = ["y", "yes"]

print()
print("Running script with the following settings:")
if answer0 not in yesList:
	print("Skipping FastQC")
elif answer0 in yesList:
	print("Using FastQC...")	
if answer1 not in yesList:
	print("Skipping mapping step...")
elif answer1 in yesList:
	print("Mapping genes with STAR aligner...")
if answer2 not in yesList:
	print("Skipping featureCounts step...")
elif answer2 in yesList:
	print("Using featureCounts step...")
print()

# Create a command log file
logFilePath = pf + "logFile.txt"
logFile = open(pf + "/" + "logFile.txt", "w")
logFile.write("Genome Used: " + g + "\n")
logFile.write("Fasta File: " + fasta_file + "\n")
logFile.write("FeatureCounts feature file: " + featuresFile + "\n\n")
logFile.write("Command Log:\n")
print("Working directory is: ", pf)

# Ensure fastq files are present. Collect sample names based on file names, & get a list of fastq files.
if answer0 in yesList or answer1 in yesList:
	fq_file_list, ff_name_actual = filecheck(rf, mf, featuresFile, gf, fasta_file, sf, ff, answer0, answer1, answer2)

# Run FastQC of fastq files
if answer0 in yesList:
	fastqc(fq_file_list, fqc, pf, ff, logFilePath)

# Run alignment of fastq files
if answer1 in yesList:
	for fq_filepair in fq_file_list:
		mapFQ(fq_filepair, ff, sg, pf, so, tf, logFilePath)

# Run featureCounts to collect read counts
if answer2 in yesList:
	featurecounts(mf, gf, rf, featuresFile)
	
# Finish
logFile.close()
print("Script complete.")
print()
exit()
##########################################################################################################