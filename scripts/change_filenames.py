# This is a tool to quickly change file names.
# This is helpful for simplifying large file names often produced as out files
# Preferred naming convention is File1_R1.fastq and File1_R1.fastq.gz

import os
print("")
print("Change filenames will alter the names of all files in the specified directory.")
print()
print("What is the target folder name/relative path?" )
print("For example, enter './' for the current directory.")
print("Alternatively, enter './fastq' if you're running this script in the parent directory of the target folder 'fastq'.")
target = input("Target path: ")
old = input("What text in the filenames should be removed/replaced? ")
new = input("What should it be replaced with? ")
print("Ok, the script will remove: '" + old + "'")
print("and replace it with: '" + new + "'")
proceed = (input("Proceed? Y/N: ")).lower()
yesList = ["y", "yes"]

if proceed in yesList:
	
	# All files
	names = (os.listdir(target))
	
	# Remove python scripts if present
	for n in names:
		if n.endswith(".py"):
			names.remove(n)
	
	newnames = []
	print("File names: ", names)
	for i in names:
		renamed = i.replace(old, new)
		print()
		print("Old name: ", i)
		print("New name: ", renamed)
		print()
		newnames.append(renamed)

	namepaths = []
	newnamepaths = []
	for n in names:
		namepaths.append(target + "/" + n)
	for n in newnames:
		newnamepaths.append(target + "/" + n)	

	ok = input("Do these all look ok? Y/N: ")
	print()
	print("Ok, renaming...")
	if ok in yesList:
		for idx in range(len(namepaths)):
			os.rename(namepaths[idx], newnamepaths[idx])

	print("Script complete.")
	print()

			
	
else:
	print("Response not valid. Exiting...")
	print()
	exit()

exit()