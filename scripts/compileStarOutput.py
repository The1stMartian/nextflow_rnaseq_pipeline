import os
####################################################
def compileStarOutput(inputFolder):
	"""
	- processes all star output files in the directory
	- From multiple star output files files, compile a single file 
	"""
	resultList = []
	fileList = os.listdir(inputFolder)

	outData = []
	rowNames = ["Sample"]

	fileNumber = 0
	for f in fileList:
		fData = []
		if f.endswith(".Log.final.out"):
			print("Extracting data from: ", f)
			fileNumber += 1
			sample = f.split(".")[0]
			fPath = inputFolder + "/" + f
			fData = [sample]
			with open(fPath, "r") as dataFile:
				for line in dataFile:
						splt = line.split("|\t")
						if fileNumber == 1:
							if len(splt) >= 2:
								#print("saving: ", splt[1].strip().strip("\n"))
								rNameDirty = splt[0].strip().strip("\n")
								rNameClean = rNameDirty.replace(",","-")	# b/c csv output format
								rowNames.append(rNameClean)
								fData.append(splt[1].strip().strip("\n"))
							else:
								continue
						elif fileNumber != 1:
							if len(splt) >= 2:
								fData.append(splt[1].strip().strip("\n"))
							else:
								continue
			# Record data
			outData.append(fData)
	return outData, rowNames

####################################################
def compileData(outFilePath, dataRows, rowNames):
	"""
	- from data and data row names, prints a single outfile (.csv)
	"""
	rowNumber = len(rowNames)
	sampleNumber = len(dataRows)

	newRows = []
	for i in range(rowNumber):
			row = rowNames[i] + ","
			for j in range(sampleNumber):
				row = row + dataRows[j][i] + ","
			newRows.append(row)

	print("Writing compiled data to: ", outFilePath)
	print()

	with open(outFilePath, "w") as out:
		for r in newRows:
			out.write(r + "\n")
	return None

####################################################
# Execute
####################################################
pwd = os.getcwd()
starDataFolder = "/mapped"
outputFileName = "compiledStarData.csv"
starDataFolderPath = pwd + starDataFolder
outputFilePath = pwd + "/" + outputFileName

compiledData, dataRowNames = compileStarOutput(starDataFolderPath)
compileData(outputFilePath, compiledData, dataRowNames)