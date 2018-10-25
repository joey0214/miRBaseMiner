#! /usr/bin/python
# python 2.7

import logging
import os
import ftplib
from Bio import SeqIO
from Bio.Seq import Seq
import editdistance
import sys
import numpy
import subprocess



def Version():
	print """miRBaseMiner Version: 0.1 """

def DMOminning():
	############################################
	## miRBaseMiner demonstration mining
	############################################
	Version()
	print "##########\n##########\nStarting Demonstration of miRBaseMiner:\n"
	# curDir = os.getcwd() ## get current work directory

	DMO_workDir = os.getcwd()  +"/DMO_miRBaseMiner"
	DMO_targetVersion = "22" ## ONLY single version
	DMO_versionList = ["20","21","22"]
	DMO_species = ["hsa"] ## can also be list ["hsa", "mmu"]
	DMO_outputFolder = "DMO_minningResults" ## if not specifiy, it will create a new folder "resultsTables" in current work dir
	DMO_NTmers = ["2","3"]  ## set it as list for furtuer extension 
	DMO_LevenshteinDistanceCutoff = 1 
	DMO_gffOptions = ["1LDmp", "ISm", "OSm"]

	print "\n===> miRBaseMiner DMO: downloading data from miRBase ftp. Version:"+", ".join(DMO_versionList)
	RetrieveMiRBase(DMO_workDir, DMO_versionList)
	
	print "\n===> miRBaseMiner DMO: Overview on annotation in miRBase releases:"+", ".join(DMO_versionList)
	Overview(DMO_workDir, DMO_versionList, DMO_outputFolder)
	
	print "\n===> miRBaseMiner DMO: Minimum free energy distribution in miRBase\nReleases:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	MFEdistribution(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species)
	
	print "\n===> miRBaseMiner DMO: Searching identical sequence in mature miRNAs.\nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	IdenticalSequenceMir(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species)
	
	print "\n===> miRBaseMiner DMO: Searching sequence overlapping in miRBase. \nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	SequenceOverlapping(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species)
	
	print "\n===> miRBaseMiner DMO: Tracking difference in miRBase releases. \nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	TrackingDiffInVersions(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species )
	
	print "\n===> miRBaseMiner DMO: Searching Reverse Complementary sequence in miRBase. \nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	ReverseComplementary(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species )
	
	print "\n===> miRBaseMiner DMO: Searching Multiple Locations for miRNAs in miRBase. \nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	MiRrMultipleLocations(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species)
	
	print "\n===> miRBaseMiner DMO: Comparsing high confidence annotation set in different miRBase release. \nVersion:"+", ".join(DMO_versionList)+"\nSpecies:"+", ".join(DMO_species)
	HighConfidenceComparison(DMO_workDir, DMO_versionList, DMO_outputFolder, DMO_species, DMO_targetVersion )
	
	print "\n===> miRBaseMiner DMO: Building sequence similarity network of miRNAs and pre-miRNAs in miRBase version:"+DMO_targetVersion+"\nfor "+", ".join(DMO_species)
	SeqSimilarityNetwork(DMO_workDir, DMO_targetVersion, DMO_outputFolder, DMO_species, DMO_LevenshteinDistanceCutoff)

	print "##########\n##########\nFinished Demonstration of miRBaseMiner\n##########\n##########\n"


## in-class function 
def makeRightPath( pathStr):
	############################################
	## check if "/" is attached to the end of the path string
	############################################
	if not pathStr.endswith("/"):
		pathStr = pathStr+"/"
	return pathStr

## public function
def RetrieveMiRBase( targetDir, versionList):
	############################################
	## downloading miRBase releases data from chosen versions into user chosen directory
	############################################
	##try to capture error of NO pigz installed
	try:
		pigzPath = subprocess.check_output("which pigz", shell=True)
		pigzPath = pigzPath.rstrip() 

		## make folder for miRBase data files
		localMiRfolder = makeRightPath(targetDir)+"miRBase"
		if not os.path.exists(localMiRfolder):
			os.makedirs(localMiRfolder)

		mirbaseFTP = "/pub/mirbase/"

		for iversion in versionList:
			print "starting downloading data of miRBase "+iversion
			iverionLocalFoldr = makeRightPath(makeRightPath(localMiRfolder)+iversion)
			## make subfolder for each miRBase release
			if not os.path.exists(iverionLocalFoldr):
				os.makedirs(iverionLocalFoldr)

			if len(os.listdir(iverionLocalFoldr)) == 0:
				## connect to ftp of miRBase
				iftpURL = mirbaseFTP +iversion+"/"
				iftp = ftplib.FTP('mirbase.org', timeout=800)
				iftp.login()
				iftp.cwd(iftpURL)
				filenames = iftp.nlst()
				## start download all files in this ftp folder
				for ifilename in filenames:
					## "database_files" and "genomes" are subfolder names
					if ifilename in {"database_files", "genomes", "database"}:
						subfolder = iverionLocalFoldr + ifilename +"/"
						iftp.cwd(iftpURL+ifilename)
						## check exist of subfolder at local, the "genomes" and "database_files" folder
						if not os.path.exists(subfolder):
							os.makedirs(subfolder)
						subFolderFiles=iftp.nlst()
						## iterate all files in subfolder
						for isubfile in subFolderFiles:
							## here try to ignore some folder that not so relevent to recent version
							if not isubfile  in {"old_releases"}:
								sublocalFilename = os.path.join(subfolder, isubfile)
								ifile = open(sublocalFilename,"wb")
								iftp.retrbinary('RETR '+ isubfile, ifile.write)
								ifile.close()
						## finished in this subfolder, now go back to previous folder
						iftp.cwd(iftpURL)
					else:
						## download files
						localFilename = os.path.join(iverionLocalFoldr,ifilename)
						ifile = open(localFilename,"wb")
						iftp.retrbinary('RETR '+ ifilename, ifile.write)
						ifile.close()
				## goodbye ftp
				iftp.quit()
				print "finished downloading data of miRBAse "+iversion
			else:
				print "NOT empty directory: skipping download data for miRBase "+iversion+"\n"
		### uncompress data
		for jversion in versionList:
			jversionFolder = makeRightPath(makeRightPath(localMiRfolder)+jversion)
			os.chdir(jversionFolder)
			## build command string
			gzCommandStr = "for file in "+jversionFolder+"*.gz; do "+pigzPath+" -d  ${file}; done"
			os.system(gzCommandStr)
		os.chdir(makeRightPath(localMiRfolder)) ## back to previous work directory
	## in case that user didn't install pigz
	except subprocess.CalledProcessError, e:
		print "Could not find executable path for pigz"
	
## in-class function 
def countingMatureInVersion( miWorkDir, miVersionList, miOutputFolder):
	############################################
	## this function is to count miRNA in each species and each version
	##  also count miRNA, summ all species in each version
	############################################
	countingSpeciesArray = [] ## in every version, how many miRNAs in each species
	sumSpeciesInVersionArray = [] ## in every version, how many miRNAs in total
	sumMatureInVersionArray = [] ## in every version, how many species in total

	## add header line 
	countingSpeciesArray.append("version\tspecies\tmatureCount\n")
	sumMatureInVersionArray.append("version\tmatureCount\n")
	sumSpeciesInVersionArray.append("version\tspeciesCount\n")
	## iterate every given versions
	for iversion in miVersionList:
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		iMatureFaSeqDict = ReadFasta2Dict(iverionFolder+"mature.fa")
		sumMatureInVersionArray.append(versionStr+"\t"+str(len(iMatureFaSeqDict))+"\n")

		speciesCountingDict = {}
		speciesList = []
		for imature in iMatureFaSeqDict:
			ispecies = imature.split("-")[0]
			speciesList.append(ispecies)
			if ispecies in speciesCountingDict:
				speciesCountingDict[ispecies] += 1
			else:
				speciesCountingDict[ispecies] = 1
		## convert dictionary into string list
		for ikey in speciesCountingDict:
			countingSpeciesArray.append(versionStr+"\t"+ ikey +"\t"+ str(speciesCountingDict[ikey])  +"\n")
		sumSpeciesInVersionArray.append(versionStr+"\t"+ str(len(list(set(speciesList)))) +"\n")
	## write into files
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	WriteList2File(countingSpeciesArray, filenamePrefix+"countingMatureInSpeciesInVersion.tsv",miOutputFolder)
	WriteList2File(sumSpeciesInVersionArray, filenamePrefix+"countingSpeciesInVersion.tsv",miOutputFolder)
	WriteList2File(sumMatureInVersionArray, filenamePrefix+"countingMatureInVersion.tsv",miOutputFolder)

## in-class function 
def matureSeqLenDistribution( miWorkDir, miVersionList, miOutputFolder):
	############################################
	## build a sequence length distribution for mature miRNAs, for each species, for each version
	############################################
	seqLengthDistArray = []
	seqLengthDistArray.append("version\tspecies\tlength\tcount\n")
	for iversion in miVersionList:
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		iMatureFaSeqDict = ReadFasta2Dict(iverionFolder+"mature.fa")
		matureLenDict = {}
		for imature in iMatureFaSeqDict:
			iSeqLength = str(len(iMatureFaSeqDict[imature]))
			ispecies = imature.split("-")[0]
			## add this ispecies into dictionary
			if not ispecies in matureLenDict:
				tmpLenDict = {}
				matureLenDict[ispecies] = tmpLenDict
			## update the length dictionary with this length
			ispeciesLenDict = matureLenDict[ispecies]
			if iSeqLength in ispeciesLenDict:
				ispeciesLenDict[iSeqLength] = ispeciesLenDict[iSeqLength] +1
			else:
				ispeciesLenDict[iSeqLength] = 1
			## updata  length dictionary of this species
			matureLenDict[ispecies] = ispeciesLenDict
		for ispecies in matureLenDict:
			iLenDict = matureLenDict[ispecies]
			for ilen in iLenDict:
				seqLengthDistArray.append(versionStr+"\t"+ispecies +"\t"+ str(ilen) +"\t"+ str(iLenDict[ilen]) +"\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	WriteList2File(seqLengthDistArray, filenamePrefix+"MatureSeqLenDistInVersions.tsv",miOutputFolder) 

## in-class function 
def buildText4SpeciesCloud( miWorkDir, miVersionList, miOutputFolder):
	############################################
	##  extract species code for build word cloud to show species abundance in mature miRNA annotation
	############################################
	for iversion in miVersionList: ## iterate version
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		iMatureFaSeqDict = ReadFasta2Dict(iverionFolder+"mature.fa")
		speciesCloudTextArray = []
		for imature in iMatureFaSeqDict:
			ispecies = imature.split("-")[0]
			speciesCloudTextArray.append(ispecies)
		## write out into file
		filenamePrefix = "miRBase_"+versionStr +"_"
		WriteList2File("  ".join(speciesCloudTextArray)+"\n",filenamePrefix+"Text4SpeciesCloud.tsv",miOutputFolder)

## in-class function 
def countingBySpecies( queryFaDict):
	############################################
	## count the amount of miRNA for each species in a given fasta dictionary
	############################################
	countingDict = {}
	for ikey in queryFaDict:
		ispecies = ikey.split("-")[0]
		if ispecies in countingDict:
			countingDict[ispecies] = countingDict[ispecies] +1
		else:
			countingDict[ispecies] = 1
	return countingDict

## in-class function 
def diffCountMatureHairpin( miWorkDir, miVersionList, miOutputFolder):
	############################################
	## this function is to count miRNA and hairpin in each species and in each version
	## to examine if there is some relationship between miRNA count and count hairpin in species or version
	############################################
	mhDiffInVersions = []
	mhDiffInVersions.append("version\tspecies\tmatureCount\thairpinCount\tdiff\n") ## header line
	mhDiffInVersionsADV = []
	mhDiffInVersionsADV.append("version\tspecies\tmatureCount\thairpinCount\tdiff\thairpinInStr\tmatureInStr\tratioInStr\n") ## header line

	for iversion in miVersionList: ## iterate version
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		iMatureFaSeqDict = ReadFasta2Dict(iverionFolder+"mature.fa") ## read in mature sequence fasta as dictionary
		iHairpinFaSeqDict = ReadFasta2Dict(iverionFolder+"hairpin.fa") ## read in hairpin sequence fasta as dictionary
		[iChildParentDict, allCHILDlist, allPARENTlist ] = ReadChildParent(iverionFolder+"miRNA.str")
		iMatureCountBySpeciesDict = countingBySpecies(iMatureFaSeqDict)
		iHairpinCountBySpeciesDict = countingBySpecies(iHairpinFaSeqDict)
		iHairpinCountBySpeciesDictCopy = iHairpinCountBySpeciesDict.copy()
		for jspecies in iMatureCountBySpeciesDict:
			countInMature = iMatureCountBySpeciesDict[jspecies]
			## in case annotated in mature without in hairpin
			if jspecies in iHairpinCountBySpeciesDict:
				countInHairpin = iHairpinCountBySpeciesDict[jspecies]
			else:
				countInHairpin = 0
			mhDiffInVersions.append(versionStr +"\t"+ jspecies +"\t"+ str(countInMature) +"\t"+ str(countInHairpin) +"\t"+ str(countInMature - countInHairpin) +"\n" )
			del iHairpinCountBySpeciesDict[jspecies]
		if len(iHairpinCountBySpeciesDict) > 0:
			for jkey in iHairpinCountBySpeciesDict:
				mhDiffInVersions.append(versionStr +"\t"+ jspecies +"\t"+ str(0) +"\t"+ str(iHairpinCountBySpeciesDict[jkey]) +"\t"+ str(0 - iHairpinCountBySpeciesDict[jkey])+"\n" )
		
		# print len(iHairpinCountBySpeciesDictCopy) 
		mergedSpecies = iMatureCountBySpeciesDict.keys()
		mergedSpecies += iHairpinCountBySpeciesDictCopy.keys()
		mergedSpecies += iChildParentDict.keys()
		print "before unique:%i, after merge %i"%(len(mergedSpecies), len(set(mergedSpecies)))
		# print iHairpinCountBySpeciesDict["hsa"]
		# print len(iHairpinCountBySpeciesDictCopy)
		# print "\t".join(iHairpinCountBySpeciesDictCopy.keys())
		missingArray = []
		for ikey in iMatureFaSeqDict.keys():
			if not ikey in allCHILDlist:
				missingArray.append(ikey)
		for jkey in iHairpinFaSeqDict.keys():
			if not jkey in allPARENTlist:
				missingArray.append(jkey)
		mergedSpecies = list(set(mergedSpecies))
		for iikey in mergedSpecies:
			imatureFaCount = 0
			ihairpinFaCount = 0
			imatureStrCount = 0
			ihairpinStrCount = 0
			ibothArmRatio = 0
			if iikey in iMatureCountBySpeciesDict.keys():
				imatureFaCount = iMatureCountBySpeciesDict[iikey]
			if iikey in iHairpinCountBySpeciesDictCopy.keys():
				ihairpinFaCount = iHairpinCountBySpeciesDictCopy[iikey]
			if iikey in iChildParentDict.keys():
				imatureStrCount = iChildParentDict[iikey][1]
				ihairpinStrCount = iChildParentDict[iikey][0]
				ibothArmRatio = float(iChildParentDict[iikey][2]) / ihairpinStrCount
			mhDiffInVersionsADV.append(versionStr +"\t"+ iikey +"\t"+ str(imatureFaCount) +"\t"+ str(ihairpinFaCount) +"\t"+ str(imatureFaCount - ihairpinFaCount) +"\t"+ str(imatureStrCount) +"\t"+ str(ihairpinStrCount) +"\t"+ str(ibothArmRatio) +"\n" )

	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	WriteList2File(mhDiffInVersions, filenamePrefix+"CountingMatureHairpinBySpeciesInVersions.tsv",miOutputFolder) 
	WriteList2File(mhDiffInVersionsADV, filenamePrefix+"CountingMatureHairpinBySpeciesInVersions-ADVANCED.tsv",miOutputFolder) 
	WriteList2File(missingArray, filenamePrefix+"CountingMatureHairpinBySpeciesInVersions-MISSING.tsv",miOutputFolder) 

## public function
def Overview( miWorkDir, miVersionList, miOutputFolder):
	############################################
	## investigate the overview of miRBase 
	## (1) the amount of miRNAs and pre-miRNAs in each species in every version
	## (2) sequence length distribution of mature miRNAs
	## (3) the varation between amount of miRNA and pre-miRNA in each species
	## (4) and word cloud for species abundancy in each version
	############################################
	## first check output folder
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	## count species and mature in versions
	countingMatureInVersion(  miWorkDir, miVersionList, miOutputFolder)
	## the distribution of mature miRNA sequence length
	matureSeqLenDistribution(miWorkDir, miVersionList, miOutputFolder)
	## extract text for building soecies cloud in each version
	buildText4SpeciesCloud(miWorkDir, miVersionList, miOutputFolder)
	# difference between the number of mature and hairpin of each species in each version
	diffCountMatureHairpin(miWorkDir, miVersionList, miOutputFolder)

	## ## ploting by R script 
	miRscriptsFolder = makeRightPath(miOutputFolder+"Rscripts")
	if not os.path.exists(miRscriptsFolder):
		os.makedirs(miRscriptsFolder)
	## ## ploting by R script --- countingMatureInVersion
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	assemblerRscript_countingMatureInVersion(filenamePrefix, miOutputFolder, miRscriptsFolder)
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"countingMatureInVersion.R")  ])
	os.system(rscriptExecuteCmd)
	## ## ploting by R script --- matureSeqLenDistribution
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	assemblerRscript_matureSeqLenDistribution(filenamePrefix, miOutputFolder, miRscriptsFolder)		
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"matureSeqLenDistribution.R")  ])
	os.system(rscriptExecuteCmd)
	## ## ploting by R script --- buildText4SpeciesCloud
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[-1])+"_"
	assemblerRscript_buildText4SpeciesCloud(filenamePrefix, miOutputFolder, miRscriptsFolder)		
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"buildText4SpeciesCloud.R")  ])
	os.system(rscriptExecuteCmd)
	## ## ploting by R script --- diffCountMatureHairpin
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"
	assemblerRscript_diffCountMatureHairpin(filenamePrefix, miOutputFolder, miRscriptsFolder)		
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"diffCountMatureHairpin.R")  ])
	os.system(rscriptExecuteCmd)

## public function
def MFEdistribution( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## investigate MFE distribution of pre-miRNA
	## also including other features related to MFE: sequence length, GC content, MFE index
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	speciesvalidatedMess = speciesValidation(miSpecies, miVersionList, miWorkDir)
	if speciesvalidatedMess[0]:
		fullDistributionMFE = []
		fullDistributionMFE.append("version\tspecies\tprecuror\tAccessionID\tMFE\thairpinLength\tGCcontent\tadjMFE\tMFEIndex\n")
		for iversion in miVersionList: ## iterate versions
			versionStr = "V"+iversion
			iverionFolder = makeRightPath(miWorkDir+iversion) 
			for ispecies in miSpecies: ## iterate species
				iHairpinFaSeqDict = ReadFasta2DictFeaturesSubset(iverionFolder+"hairpin.fa", ispecies) ## read fasta sequence
				iMFEhairpinDict = ReadMFEfromMiRStr(iverionFolder+"miRNA.str") ## read in MFE profile
				for ihairpin in iHairpinFaSeqDict: ## iterate each pre-miRNA
					if ihairpin in iMFEhairpinDict: ## check if this pre-miRNA has MFE profile
						iAccessionID = iHairpinFaSeqDict[ihairpin].split("\t")[0] ##  the accession number of pre-miRNA
						iHairpinLength = float(iHairpinFaSeqDict[ihairpin].split("\t")[1]) ## pre-miRNA sequence length
						iHairpinGC = float(iHairpinFaSeqDict[ihairpin].split("\t")[2]) ## GC content
						iHairpinMFE = float(iMFEhairpinDict[ihairpin]) ## MFE
						iAdjMFE = iHairpinMFE / iHairpinLength ## adjusted MFE
						iMFEindex = iAdjMFE / iHairpinGC ## MFE index 
						## add to output list 
						fullDistributionMFE.append(versionStr +"\t"+ ispecies +"\t"+ ihairpin +"\t"+ iAccessionID +"\t"+ str(iHairpinMFE) +"\t"+ str(iHairpinFaSeqDict[ihairpin].split("\t")[1]) +"\t"+ str(iHairpinGC) +"\t"+ str(iAdjMFE) +"\t"+ str(iMFEindex) +"\n")
					else:
						print "Couldn't find MFE for %s in miRNA.str!"%(ihairpin )
		## write out into file
		filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
		WriteList2File(fullDistributionMFE, filenamePrefix + "fullDistribution4MFE.tsv",miOutputFolder)

		## ## ploting by R script --- MFEdistribution
		miRscriptsFolder = makeRightPath(miOutputFolder+"Rscripts")
		if not os.path.exists(miRscriptsFolder):
			os.makedirs(miRscriptsFolder)
		filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
		assemblerRscript_MFEdistribution(filenamePrefix, miOutputFolder, miRscriptsFolder)
		rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"MFEdistribution.R")  ])
		os.system(rscriptExecuteCmd)
	else:
		## species validation failed
		print speciesvalidatedMess[1]
	
## in-class function 
def identicalWithinSpecies( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## find miRNAs share identical sequence within each given species
	############################################
	uniqueSeqWithINSpeciesArray = []
	uniqueSeqWithINSpeciesArray.append("version\tspecies\tmolecular\tsequence\tcount\tnameIDs\n") ## add header line
	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for ispecies in miSpecies: ## iterate species 
			for imolecular in ["mature", "hairpin"]: ## iterate molecular types, automatically check mature and hairpin sequence
				iUniqueSequenceDict = {}
				iFastaSeqDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+imolecular+".fa", ispecies)
				for imatureID in iFastaSeqDict:
					imatureSeq = iFastaSeqDict[imatureID]
					if imatureSeq in iUniqueSequenceDict:
						iUniqueSequenceDict[imatureSeq] = iUniqueSequenceDict[imatureSeq] +";"+imatureID
					else:
						iUniqueSequenceDict[imatureSeq] = imatureID
				## convert dictionary into array
				for ikey in iUniqueSequenceDict:
					if ";" in iUniqueSequenceDict[ikey]:
						uniqueSeqWithINSpeciesArray.append(versionStr +"\t"+ ispecies +"\t"+ imolecular +"\t"+ ikey +"\t"+ str(iUniqueSequenceDict[ikey].count(";")+1)+"\t"+ iUniqueSequenceDict[ikey] +"\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(uniqueSeqWithINSpeciesArray, filenamePrefix+"IdenticalSequenceWithInSpecies.tsv",miOutputFolder) 

	## ## ploting by R script --- MFEdistribution
	miRscriptsFolder = makeRightPath(miOutputFolder+"Rscripts")
	if not os.path.exists(miRscriptsFolder):
		os.makedirs(miRscriptsFolder)
	assemblerRscript_identicalWithinSpecies(filenamePrefix, miOutputFolder, miRscriptsFolder)
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"identicalWithinSpecies.R")  ])
	os.system(rscriptExecuteCmd)

## in-class function 
def identicalBetweenSpecies( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## find miRNAs share identical sequence between given species 
	############################################
	uniqueSeqBetweenSpeciesArray = []
	uniqueSeqBetweenSpeciesArray.append("version\tmolecular\tsequence\t"+"\t".join(miSpecies)+ "\n") ## add header line
	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for imolecular in ["mature", "hairpin"]: ## iterate molecular types, automatically check mature and hairpin sequence
			uniqueSeqBetweenSpecies = {}
			for i in range(len(miSpecies)):
				ispecies = miSpecies[i]
				iFastaSeqDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+imolecular+".fa", ispecies)
				for ifastaID in iFastaSeqDict:
					ifastaSeq = iFastaSeqDict[ifastaID]
					if ifastaSeq in uniqueSeqBetweenSpecies:
						tmparr = uniqueSeqBetweenSpecies[ifastaSeq]
						tmparr[i].append(ifastaID) 
						uniqueSeqBetweenSpecies[ifastaSeq] = tmparr
					else:
						tmparr = [[] for _ in range(len(miSpecies))]
						tmparr[i].append(ifastaID)
						uniqueSeqBetweenSpecies[ifastaSeq] = tmparr
			for iseqkey in uniqueSeqBetweenSpecies:
				tmplist = []
				for j in range(len(miSpecies)):
					tmplist.append(";".join(uniqueSeqBetweenSpecies[iseqkey][j]))
				uniqueSeqBetweenSpeciesArray.append(versionStr +"\t"+  imolecular +"\t"+ iseqkey +"\t"+  "\t".join(tmplist) +"\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(uniqueSeqBetweenSpeciesArray, filenamePrefix+"IdenticalSequenceBetweenSpecies.tsv",miOutputFolder) 

## public function
def IdenticalSequenceMir( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## investigate full length identical sequence in miRNA/pre-miRNA, share sequence between different species
	## require more than one species in miSpecies
	############################################
	## prepare output folder and work directory
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	identicalWithinSpecies(miWorkDir, miVersionList, miOutputFolder, miSpecies)
	## check if miSpecies fulfill the requirement
	if len(miSpecies) > 1 :
		identicalBetweenSpecies(miWorkDir, miVersionList, miOutputFolder, miSpecies)

## in-class function 
def buildName2NameAccessiondict( fastaDict1, fastaDict2):
	############################################
	## merge keys from two fasta dictionaries; 
	## build a new dictionary, ID --- ID/accession
	############################################
	name2nameAccession = {}
	for ikey in fastaDict1.keys():
		name2nameAccession[ikey.split("/")[0]] = ikey
	for jkey in fastaDict2.keys():
		name2nameAccession[jkey.split("/")[0]] = jkey
	return name2nameAccession

## in-class function
def findingSequenceOverlapping( iFastaSeqDict, versionStr, ispecies, imolecular,mirParentChildDict, name2nameAccessionDict, hairpinFaDict):
	############################################
	## find sequence overlapping with others in a given fasta dictionary
	############################################
	overlapedPairs = []
	matureOverlappingArray = []
	hairpinOverlappingArray = []
	matureKeylist = iFastaSeqDict.keys()
	for i in range(0, len(matureKeylist) ):
		ifastaSeq = matureKeylist[i]
		for j in range(i+1, len(matureKeylist) ):
			jfastaSeq = matureKeylist[j]
			ISoverlap = False
			if iFastaSeqDict[ifastaSeq] in  iFastaSeqDict[jfastaSeq] :
				ISoverlap = True
			elif iFastaSeqDict[jfastaSeq] in  iFastaSeqDict[ifastaSeq]:
				ISoverlap = True
			if ISoverlap:
				if len(overlapedPairs) == 0:
					overlapedPairs.append([ifastaSeq,jfastaSeq])
				else:
					ISIncluded = False
					includedIndex = -1
					for ipairs in overlapedPairs:
						## check if match with any exist pairs
						if ifastaSeq in ipairs:
							includedIndex = overlapedPairs.index(ipairs)
							ISIncluded = True
						elif jfastaSeq in ipairs:
							includedIndex = overlapedPairs.index(ipairs)
							ISIncluded = True
					## add as new pair
					if ISIncluded:
						includedPair = overlapedPairs[includedIndex]
						overlapedPairs.pop(includedIndex)
						includedPair.append(ifastaSeq)
						includedPair.append(jfastaSeq)
						newipairs = list(set(includedPair))
						overlapedPairs.append(newipairs)
					else:
						overlapedPairs.append([ifastaSeq,jfastaSeq])
	pairIndex = 1
	for collectedPairs in overlapedPairs:
		for i in collectedPairs:
			faSeqStr = iFastaSeqDict[i]
			parentChildNames = mirParentChildDict[i.split("/")[0]]
			if imolecular == "mature":
				if ";" in  parentChildNames:
					for iname in parentChildNames.split(";"):
						iname = name2nameAccessionDict[iname]
						matureOverlappingArray.append(versionStr+"\t"+ ispecies+"\t"+ i +"\t"+ faSeqStr+"\t"+iname+"\t"+hairpinFaDict[iname]+"\n")
				else:
					parentChildNames = name2nameAccessionDict[parentChildNames]
					matureOverlappingArray.append(versionStr+"\t"+ ispecies+"\t"+ i +"\t"+ faSeqStr+"\t"+parentChildNames+"\t"+hairpinFaDict[parentChildNames]+"\n")
			else:
				if ";" in  parentChildNames:
					tmpNameList = []
					for iname in parentChildNames.split(";"):
						iname = name2nameAccessionDict[iname]
						tmpNameList.append(iname)
					hairpinOverlappingArray.append(versionStr+"\t"+ ispecies+"\t"+ i+"\t"+ faSeqStr+"\t"+";".join(tmpNameList)+"\n")
				else:
					parentChildNames = name2nameAccessionDict[parentChildNames]
					hairpinOverlappingArray.append(versionStr+"\t"+ ispecies+"\t"+ i+"\t"+ faSeqStr+"\t"+parentChildNames+"\n")
		if imolecular == "mature":	
			matureOverlappingArray.append(" \n")
		else:
			hairpinOverlappingArray.append(" \n")
		pairIndex += 1
	if imolecular == "mature":	
		return matureOverlappingArray
	else:
		return hairpinOverlappingArray

## public function
def SequenceOverlapping( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## find sequence overlapping with others in given species list and version list
	############################################
	## prepare output folder and work directory
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	matureOverlappingArray = []
	matureOverlappingArray.append("version\tspecies\tnameAccession\tsequence\tprecursor\thairpinSequence\n\n")
	hairpinOverlappingArray = []
	hairpinOverlappingArray.append("version\tspecies\tnameAccession\tsequence\tmatureAccession\n\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		mirParentChildDict = ReadParentChildFromStr(iverionFolder+"miRNA.str", miSpecies)
		for ispecies in miSpecies: ## iterate species 
			## since some pre-miRNAs or miRNAs DO NOT have genome coordinates in gff file
			## so here we extract the relationship between miRNA and pre-miRNA from the file "miRNA.str"
			hairpinFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"hairpin.fa", ispecies)
			matureFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"mature.fa", ispecies)
			## ID --- ID/Accession dictionary
			name2nameAccessionDict = buildName2NameAccessiondict(hairpinFaDict, matureFaDict )
			for imolecular in ["mature", "hairpin"]: ## iterate molecular types, automatically check mature and hairpin sequence
				# iUniqueSequenceDict = {}
				# iFastaSeqDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+imolecular+".fa", ispecies)
				if imolecular == "mature":
					iFastaSeqDict = matureFaDict
				else:
					iFastaSeqDict = hairpinFaDict
				iOverlappingArray = findingSequenceOverlapping(iFastaSeqDict, versionStr, ispecies, imolecular,mirParentChildDict, name2nameAccessionDict, hairpinFaDict)
				if imolecular == "mature":
					matureOverlappingArray = matureOverlappingArray + iOverlappingArray
				else:
					hairpinOverlappingArray = hairpinOverlappingArray+ iOverlappingArray
	## write file 
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(matureOverlappingArray, filenamePrefix+"SequenceOverlapping_mature.tsv",miOutputFolder) 
	WriteList2File(hairpinOverlappingArray, filenamePrefix+"SequenceOverlapping_hairpin.tsv",miOutputFolder) 

## in-class function 
def endNucleotide( miWorkDir, miVersionList, miOutputFolder, miSpecies, miNTmers):
	############################################
	## investigate the nucleotide combination at both termini of miRNAs
	## the nucleotide number is determined by miNTmers
	############################################
	endNucleotideCombinationArray = []
	endNucleotideCombinationArray.append("version\tspecies\tendPosition\tNucleotides\tcount\tratio\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		# mirParentChildDict = ReadParentChildFromStr(iverionFolder+"miRNA.str", miSpecies)
		for ispecies in miSpecies: ## iterate species 
			iMatureFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"mature.fa", ispecies)
			ntMers5pDict = {}
			ntMers3pDict = {}
			for iMers in  miNTmers:
				iMers = int(iMers)
				for ifastaID in  iMatureFaDict:
					ifastaSeq = iMatureFaDict[ifastaID]
					iNT5mer = ifastaSeq[: iMers]
					iNT3mer = ifastaSeq[-iMers: ]
					## 5 prime end
					if iNT5mer in ntMers5pDict.keys():
						ntMers5pDict[iNT5mer] += 1
					else:
						ntMers5pDict[iNT5mer] = 1
					## 3 prime end 
					if iNT3mer in ntMers3pDict.keys():
						ntMers3pDict[iNT3mer] += 1
					else:
						ntMers3pDict[iNT3mer] = 1 
			## calculate ration for each nucleotide combination
			for ikey in ntMers5pDict.keys():
				iRatio = float(ntMers5pDict[ikey]) / len(iMatureFaDict)
				endNucleotideCombinationArray.append(versionStr +"\t"+ ispecies +"\t"+"5p"+str(len(ikey))+"\t"+ ikey +"\t"+ str(ntMers5pDict[ikey]) +"\t"+ str(iRatio) +"\n" )
			for jkey in ntMers3pDict.keys():
				jRatio = float(ntMers3pDict[jkey]) / len(iMatureFaDict)
				endNucleotideCombinationArray.append(versionStr +"\t"+ ispecies +"\t"+"3p"+str(len(jkey))+"\t"+ jkey +"\t"+ str(ntMers3pDict[jkey]) +"\t"+ str(jRatio) +"\n" )
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(endNucleotideCombinationArray, filenamePrefix+"EndNucleotideCombination_mature.tsv",miOutputFolder) 

## in-class function 
def polyAtailing( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## investigate the poly(A) at 3' end of miRNA
	## ONLY 3' end
	## ONLY poly(A)
	############################################
	polyAtailingArray = []
	polyAtailingArray.append("version\tspecies\ttype\tmiRNA_ID\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for ispecies in miSpecies: ## iterate species 
			iMatureFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"mature.fa", ispecies)

			for ifastaID in iMatureFaDict.keys():
				ifastaSeq = iMatureFaDict[ifastaID] 
				for j in range(1,6): ## set 5 as max of polyA tailing, iterate each possibility
					for theNT in ["A"]: ## for later extension
						ntPattern = theNT * j
						if ifastaSeq.endswith(ntPattern):
							polyAtailingArray.append(versionStr +"\t"+ ispecies +"\t"+ str(j)+theNT +"\t"+ ifastaID + "\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(polyAtailingArray, filenamePrefix+"polyAtailing_mature.tsv",miOutputFolder)

## in-class function 
def polyNTheading( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## investigate the polynucleotide at 5' end of miRNA
	## ONLY 5' end
	## add code for supporting polynucleotide at 3' end of miRNA
	############################################
	polyNTheadingArray = []
	polyNTheadingArray.append("version\tspecies\tterminal\tpolyNucleotide\tmiRNA_ID\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for ispecies in miSpecies: ## iterate species 
			iMatureFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"mature.fa", ispecies)
			for ifastaID in iMatureFaDict.keys():
				ifastaSeq = iMatureFaDict[ifastaID] 
				for j in range(1,6): ## set 5 as max of polyA tailing
					for theNT in ["A","G","C","U","T"]: ## for later extension
						ntPattern = theNT * j
						if ifastaSeq.startswith(ntPattern):
							polyNTheadingArray.append(versionStr +"\t"+ ispecies +"\t5p\t"+  str(j)+theNT +"\t"+ ifastaID + "\n")
						# elif ifastaSeq.endswith(ntPattern):
						# 	polyNTheadingArray.append(versionStr +"\t"+ ispecies +"\t3p\t"+  str(j)+theNT +"\t"+ ifastaID + "\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	WriteList2File(polyNTheadingArray, filenamePrefix+"polyNTheading_mature.tsv",miOutputFolder)

## in-class function
def seedRegionSeq( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## build seed region sequence profile
	## for (1) SeqLogo; (2) heat map of similarity in seed region
	############################################
	seedRegionSeqAlign = []
	seedRegionSeqAlign.append("version\tspecies\ttype\tmiRNA_ID\n")
	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for ispecies in miSpecies: ## iterate species 
			iMatureFaDict = ReadFasta2DictSubsetKeepUFullKey(iverionFolder+"mature.fa", ispecies)
			seedSeqArray = []
			for ifastaID in iMatureFaDict:
				ifastaSeq = iMatureFaDict[ifastaID]
				seedSeqArray.append(ifastaID +"\t"+ ifastaSeq[1:7] +"\n") ## canonical seed region sequence
			## write out into file
			filenamePrefix = "miRBase_"+versionStr +"_"+ispecies+"_"
			WriteList2File(seedSeqArray, filenamePrefix+"SeedRedionSequence.tsv",miOutputFolder)

## public function
def NucleotideCombination( miWorkDir, miVersionList, miOutputFolder, miSpecies, miNTmers):
	############################################
	## investigate the nucleotide combination of miRNA
	## ONLY available for mature miRNAs
	## (1) endNucleotide; (2) polyAtailing; (3) polyNTheading; (4) seedRegionSeq.
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	# ## part X nts at both 5 and 3 prime end
	endNucleotide(miWorkDir, miVersionList, miOutputFolder, miSpecies, miNTmers)
	## check polyA tail
	polyAtailing(miWorkDir, miVersionList, miOutputFolder, miSpecies)
	## check polyNT header
	polyNTheading(miWorkDir, miVersionList, miOutputFolder, miSpecies)
	## canonical seed region 
	seedRegionSeq(miWorkDir, miVersionList, miOutputFolder, miSpecies)

	## ## ploting by R script --- MFEdistribution
	miRscriptsFolder = makeRightPath(miOutputFolder+"Rscripts")
	if not os.path.exists(miRscriptsFolder):
		os.makedirs(miRscriptsFolder)
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])+"_"+"-".join(miSpecies)+"_"
	assemblerRscript_NucleotideCombination(filenamePrefix, miOutputFolder, miRscriptsFolder)
	rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+"NucleotideCombination.R")  ])
	os.system(rscriptExecuteCmd)

## 	public function
def TrackingDiffInVersions( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## this function is to trick what happen in "miRNA.diff" file
	## this function can catch some entries that have multiple change tags (i.e. NAME and SEQUENCE)  
	## summarize :
	## ## (1)how many NEW miRNAs were added  
	## ## (2)how many miRNAs were  DELETE 
	## ## (3)how many miRNAs were  changed SEQUENCE
	## ## (4)how many miRNAs were  changed NAME

	## ## # NEW           new to the database
	## ## # DELETE        removed since the last release
	## ## # SEQUENCE      sequence has changed
	## ## # NAME          ID has changed
	############################################
	## prepare output folder and work directory
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	
	subsetDiffArray = []
	subsetDiffArray.append("version\tspecies\taccession\tmolecularName\tchangeType\n")
	countingDiffArray = []
	countingDiffArray.append("version\tmolecularType\tspecies\tchangeType\tcount\n")
	sequenceChangeArray = []
	sequenceChangeArray.append("version\tspecies\tmolecularType\tnameAccession\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		iDiffFile = iverionFolder+"miRNA.diff"
		diffHandle = open(iDiffFile,"r")
		matureCountingDict = {}
		hairpinCountingDict = {}
		for line in diffHandle:
			if not line.startswith("#"):
				line = line.rstrip()
				linesplit = line.split("\t")
				keyStr = linesplit[-1]
				accessionStr = linesplit[0]
				nameIDstr = line.split("\t")[1]
				ispecies =  nameIDstr.split("-")[0]
				### count how many miRNAs/pre-miRNAs had changed
				if ispecies in miSpecies:
					if accessionStr.startswith("MIMAT"):
						## this mean the entry line is talking about miRNA
						if len(linesplit) >3:
							print "%s has multiple change in version %s"%(nameIDstr, versionStr)
							for ikey in linesplit[2:]:
								newkey = ispecies+"-"+ikey
								if newkey in matureCountingDict:
									matureCountingDict[newkey] += 1
								else:
									matureCountingDict[newkey] = 1
								if ikey == "SEQUENCE":
									sequenceChangeArray.append("\t".join([versionStr, ispecies, "mature", (nameIDstr+"/"+accessionStr)]) +"\n")
						else:
							newkey = ispecies+"-"+keyStr
							if newkey in matureCountingDict:
								matureCountingDict[newkey] += 1
							else:
								matureCountingDict[newkey] = 1
							if keyStr == "SEQUENCE":
								sequenceChangeArray.append("\t".join([versionStr, ispecies, "mature", (nameIDstr+"/"+accessionStr)])+"\n")
					else:
						## this mean the entry line is talking about pre-miRNA
						if len(linesplit) >3:
							print "%s has multiple change in version %s"%(nameIDstr, versionStr)
							for ikey in linesplit[2:]:
								newkey = ispecies+"-"+ikey
								if newkey in matureCountingDict:
									matureCountingDict[newkey] += 1
								else:
									matureCountingDict[newkey] = 1 
								if ikey == "SEQUENCE":
									sequenceChangeArray.append("\t".join([versionStr, ispecies, "hairpin", (nameIDstr+"/"+accessionStr)])+"\n")
						else:
							newkey = ispecies+"-"+keyStr
							if newkey in hairpinCountingDict:
								hairpinCountingDict[newkey] += 1
							else:
								hairpinCountingDict[newkey] = 1
							if keyStr == "SEQUENCE":
								sequenceChangeArray.append("\t".join([versionStr, ispecies, "hairpin", (nameIDstr+"/"+accessionStr)])+"\n")
					##
					subsetDiffArray.append(versionStr +"\t"+ ispecies +"\t"+ line+"\n")
		diffHandle.close()
		for iKey in matureCountingDict:
			countingDiffArray.append(versionStr +"\t"+ "mature" +"\t"+ "\t".join(iKey.split("-")) +"\t"+ str(matureCountingDict[iKey]) +"\n")
		for jKey in hairpinCountingDict:
			countingDiffArray.append(versionStr +"\t"+ "hairpin" +"\t"+ "\t".join(jKey.split("-")) +"\t"+ str(hairpinCountingDict[jKey]) +"\n" )
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1])  +"_"+"-".join(miSpecies)+"_"
	WriteList2File(countingDiffArray, filenamePrefix+"CountingDiffInVersions.tsv",miOutputFolder)
	WriteList2File(subsetDiffArray , filenamePrefix+"DiffTableInVersions.tsv",miOutputFolder)
	WriteList2File(sequenceChangeArray , filenamePrefix+"SEQUENCEchangeInVersions.tsv",miOutputFolder)

## in-class function
def findingReverseComplementary( queryFasta,versionStr,ispecies, iMolecular):
	############################################
	## finding Reverse Complementary in a given fasta dictionary
	############################################
	## fold into unqie sequence dictionary 
	iUniqueFaDict = {}
	reverseComplementaryMiRArray = []
	fullLenRClist = []
	overlappingRClist = []

	for ifastaID in queryFasta:
		ifastaSeq = queryFasta[ifastaID].replace("U", "T")
		if ifastaSeq in iUniqueFaDict.keys():
			iUniqueFaDict[ifastaSeq] = iUniqueFaDict[ifastaSeq] + ";"+ ifastaID
		else:
			iUniqueFaDict[ifastaSeq] =  ifastaID
	
	querySeqStrList = iUniqueFaDict.keys()
	for i in range(0, len(querySeqStrList)):
		iQuerySeq = querySeqStrList[i]
		iQueryAccId = iUniqueFaDict[iQuerySeq]
		iQuerySeqRC = Seq(iQuerySeq)
		iQuerySeqRC = iQuerySeqRC.reverse_complement()

		for j in range(i+1, len(querySeqStrList)):
			jQuerySeq = querySeqStrList[j]
			jQueryAccID = iUniqueFaDict[jQuerySeq]

			if str(iQuerySeqRC) == jQuerySeq:
				reverseComplementaryMiRArray.append(versionStr +"\t"+ ispecies +"\t"+ iMolecular +"\t"+ iQueryAccId +"\t"+ iQuerySeq +"\t"+ jQueryAccID +"\t"+ jQuerySeq +"\t"+ "FULL" +"\n")
				fullLenRClist.append(iQueryAccId)
				fullLenRClist.append(jQueryAccID)
			elif len(str(iQuerySeqRC)) > len(jQuerySeq):
				if jQuerySeq in str(iQuerySeqRC):
					reverseComplementaryMiRArray.append(versionStr +"\t"+ ispecies +"\t"+ iMolecular +"\t"+ iQueryAccId +"\t"+ iQuerySeq +"\t"+ jQueryAccID +"\t"+ jQuerySeq +"\t"+ "OVERLAP" +"\n")
					overlappingRClist.append(iQueryAccId)
					overlappingRClist.append(jQueryAccID)
			elif len(str(iQuerySeqRC)) < len(jQuerySeq):
				if str(iQuerySeqRC) in jQuerySeq:
					reverseComplementaryMiRArray.append(versionStr +"\t"+ ispecies +"\t"+ iMolecular +"\t"+ iQueryAccId +"\t"+ iQuerySeq +"\t"+ jQueryAccID +"\t"+ jQuerySeq +"\t"+ "OVERLAP" +"\n")
					overlappingRClist.append(iQueryAccId)
					overlappingRClist.append(jQueryAccID)
	return [list(set(fullLenRClist)) ,list(set(overlappingRClist)),  reverseComplementaryMiRArray]

## public function 
def ReverseComplementary( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## reverse complementary within each given species in given versions
	## only focus on the sequence, regardless the sequence annotated on which strand
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	
	reverseComplementaryMiRArray = []
	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for iMolecular in ["mature", "hairpin"]: ## iterate molecular types
			for ispecies in miSpecies: ## iterate species
				iFastaSeqDict = ReadFasta2DictSubsetKeepUFullKey((iverionFolder+iMolecular+".fa"), ispecies)
				[fullLenRClist ,overlappingRClist,  iReverseComplementaryMiRArray] = findingReverseComplementary( iFastaSeqDict,versionStr,ispecies, iMolecular)
				reverseComplementaryMiRArray += iReverseComplementaryMiRArray
	## write file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1]) +"_"+"-".join(miSpecies)+"_"
	WriteList2File(reverseComplementaryMiRArray, filenamePrefix+"ReverseComplementaryInVersions.tsv",miOutputFolder)
	
## public function
def MiRrMultipleLocations( miWorkDir, miVersionList, miOutputFolder, miSpecies ):
	############################################
	## collect genome coordinates for miRNA/pre-miRNA from gff annotation file
	## use miRNA/pre-miRNA ID/Accession as key to pull multiple location information together
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	
	multiLocMiRArrayInfoArray = []
	multiLocMiRArrayInfoArray.append("version\tspecies\tmolecularType\tmolecularName\tlocationCount\tlocations\n")

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 
		for ispecies in miSpecies: ## iterate species
			iGFFpath= iverionFolder +"genomes/"+ ispecies +".gff3"
			## if gff3 is not available in this version, check if gff2 is available
			if not os.path.isfile(iGFFpath):
				iGFFpath = iGFFpath.replace("gff3", "gff2")
				if not os.path.isfile(iGFFpath):
					iGFFpath = iGFFpath.replace("gff2", "gff")
			## if still could not find gff file, maybe in this version, the gff file is not available for this species
			if not os.path.isfile(iGFFpath):
				sys.exit("Could not find gff file!")
				print "Could not find %s! Check miRBase ftp site, if gff file is available for this species %s"%(iGFFpath, ispecies)
			## seperate gff parser gor gff3 and gff2/gff
			if iGFFpath.endswith("gff3"):
				[matureEntryCount, hairpinEntryCount, matureMultiLocationDict, hairpinMultiLocationDict] = ReadGffV3GetGenomeCoordinates(iGFFpath)
			else:
				[matureEntryCount, hairpinEntryCount, matureMultiLocationDict, hairpinMultiLocationDict] = ReadGffV2GetGenomeCoordinates(iGFFpath)

			print "%s miRNA entries are annotated, %s hairpin entries are annotated in %s" %(str(matureEntryCount), str(hairpinEntryCount), iversion)
			## convert location dictionary into list for mature miRNA
			for ikey in matureMultiLocationDict:
				multiLocMiRArrayInfoArray.append(versionStr +"\t"+ ispecies +"\t"+ "mature" +"\t"+ ikey+"\t"+ str(matureMultiLocationDict[ikey].count(";")+1) +"\t"+matureMultiLocationDict[ikey] +"\n")
			## convert location dictionary into list for mature pre-miRNA
			for jkey in hairpinMultiLocationDict:
				multiLocMiRArrayInfoArray.append(versionStr +"\t"+ ispecies +"\t"+ "hairpin" +"\t"+ jkey+"\t"+ str(hairpinMultiLocationDict[jkey].count(";")+1) +"\t"+hairpinMultiLocationDict[jkey] +"\n")
	## write out into file
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1]) +"_"+"-".join(miSpecies)+"_"
	WriteList2File(multiLocMiRArrayInfoArray, filenamePrefix+"MultiLocationsInVersions.tsv",miOutputFolder)

## public function
def HighConfidenceComparison( miWorkDir, miVersionList, miOutputFolder, miSpecies, miQueryVersion):
	############################################
	## investigate high confidence files in version 20, 21 an 22
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	## read in from dat file
	highConfV20 = ReadAC2IDfromDat_Full(( makeRightPath(miWorkDir+"20")+ "miRNA_high_conf.dat"))
	highConfV21 = ReadAC2IDfromDat_Full(( makeRightPath(miWorkDir+"21")+ "miRNA_high_conf.dat"))
	highConfV22 = ReadAC2IDfromDat_Full(( makeRightPath(miWorkDir+"22")+ "miRNA_high_conf.dat"))
	
	## build ID/Accession list, joined from three versions
	precursorUniqueAccessionList = list(set(highConfV20[0] + highConfV21[0] + highConfV22[0]))
	matureUniqueAccessionList    = list(set(highConfV20[2] + highConfV21[2] + highConfV22[2]))

	## pull them together to build a table 
	precursorHighConfProfileArray = []
	precursorHighConfProfileArray.append("Accession\tmiRBase_V20\tmiRBase_V21\tmiRBase_V22\n")
	for iaccession in precursorUniqueAccessionList:
		presenceInV20 = "NONE"
		presenceInV21 = "NONE"
		presenceInV22 = "NONE"
		## check labeled with "high_conf" in which version
		if iaccession in highConfV20[1].keys(): ## check in version 20
			presenceInV20 =  "highconf" # highConfV20[1][iaccession]
		if iaccession in highConfV21[1].keys(): ## check in version 21
			presenceInV21 =  "highconf" # highConfV21[1][iaccession]
		if iaccession in highConfV22[1].keys(): ## check in version 22
			presenceInV22 =  "highconf" # highConfV22[1][iaccession]
		## 
		if presenceInV20 == "highconf":
			if presenceInV21 == "highconf":
				if presenceInV22 != "highconf":
					presenceInV22 = "REMOVED"
			else:
				presenceInV21 =  "REMOVED"
				if presenceInV22 != "highconf":
					presenceInV22 = "REMOVED"
		else:
			if presenceInV21 == "highconf":
				if presenceInV22 != "highconf":
					presenceInV22 = "REMOVED"

		precursorHighConfProfileArray.append(iaccession +"\t"+ presenceInV20 +"\t"+ presenceInV21 +"\t"+ presenceInV22 +"\n")
	## write out into file
	filenamePrefix = "miRBase_v20-v21-v22"+"_"
	WriteList2File(precursorHighConfProfileArray, filenamePrefix+"HighConfidenceProfile.tsv", miOutputFolder)

	for ispecies in miSpecies:
		speciesArray = []
		for ientry in precursorHighConfProfileArray:
			if (ispecies+"-") in ientry:
				speciesArray.append(ientry)
		## write file for each species
		filenamePrefix = "miRBase_v20-v21-v22"+"_"+ispecies+"_"
		WriteList2File(speciesArray, filenamePrefix+"HighConfidenceProfile.tsv", miOutputFolder)
	
	##### =====> mature miRNA 
	matureHighConfProfileArray = []
	matureHighConfProfileArray.append("Accession\tmiRBase_V20\tmiRBase_V21\tmiRBase_V22\n")
	for jaccession in matureUniqueAccessionList:
		presenceInV20m = "NONE"
		presenceInV21m = "NONE"
		presenceInV22m = "NONE"
		## check labeled with "high_conf" in which version
		if jaccession in highConfV20[3].keys(): ## check in version 20
			presenceInV20m = "highconf" # highConfV20[3][jaccession]
		if jaccession in highConfV21[3].keys(): ## check in version 21
			presenceInV21m = "highconf" # highConfV21[3][jaccession]
		if jaccession in highConfV22[3].keys(): ## check in version 22
			presenceInV22m = "highconf" # highConfV22[3][jaccession]
		## 
		if presenceInV20m == "highconf":
			if presenceInV21m == "highconf":
				if presenceInV22m != "highconf":
					presenceInV22m = "REMOVED"
			else:
				presenceInV21m =  "REMOVED"
				if presenceInV22m != "highconf":
					presenceInV22m = "REMOVED"
		else:
			if presenceInV21m == "highconf":
				if presenceInV22m != "highconf":
					presenceInV22m = "REMOVED"
	
		matureHighConfProfileArray.append(jaccession +"\t"+ presenceInV20m +"\t"+ presenceInV21m +"\t"+ presenceInV22m +"\n")
	## write out into file
	filenamePrefix = "miRBase_v20-v21-v22"+"_"
	WriteList2File(matureHighConfProfileArray, filenamePrefix+"mature.HighConfidenceProfile.tsv", miOutputFolder)
	mergedMatureDict = (highConfV20[3]).copy()
	mergedMatureDict.update(highConfV21[3])
	mergedMatureDict.update(highConfV22[3])
	for ispecies in miSpecies:
		speciesArray = []
		for ientry in matureHighConfProfileArray:
			if "Accession" in ientry:
				speciesArray.append(ientry)
			else:
				if (ispecies.lower()+"-") in mergedMatureDict[ ientry.split("\t")[0] ]:
					speciesArray.append(ientry)
		## write file for each species
		filenamePrefix = "miRBase_v20-v21-v22"+"_"+ispecies+"_"
		WriteList2File(speciesArray, filenamePrefix+"mature.HighConfidenceProfile.tsv", miOutputFolder)


## in-class function
def calculatePairwiseLevenshteinDistance( iFastaSeqDict, imolecular, miLevDist):
	############################################
	## based on a given fasta dictionary, calculate the Levenshtein Distance of each two miRNAs
	############################################
	name2NameAccessionDict = {}
	distArray = [] ## full distance table
	minimumDistArray = [] ## subset distance table
	subsetLevDistList = [] ## only ID/Accession list
	
	minimumDistDict = {}
	for ikey in iFastaSeqDict.keys():
		minimumDistDict[ikey] = 999
	fastaIDAccList = iFastaSeqDict.keys()
	for i in range(0,len(fastaIDAccList)):
		ifastaID = fastaIDAccList[i]
		ifastaSeq = iFastaSeqDict[ifastaID]
		name2NameAccessionDict[ifastaID.split("/")[0]] = ifastaID
		mindist = 999
		for j in range((i+1),len(fastaIDAccList)):
			jfastaID = fastaIDAccList[j]
			jfastaSeq = iFastaSeqDict[jfastaID]
			idist = editdistance.eval(ifastaSeq, jfastaSeq)
			## if this distance within the  Levenshtein Distance cutoff, then add into output list 
			if idist <= miLevDist: 
				distArray.append(ifastaID +"\t"+ jfastaID +"\t"+ str(idist) +"\n")
			## checl if this disyance is the minimum distance of this miRNA or pre-miRNA
			if idist < minimumDistDict[ifastaID]:
				minimumDistDict[ifastaID] = idist
			if idist < minimumDistDict[jfastaID]:
				minimumDistDict[jfastaID] = idist
	for ikey in minimumDistDict.keys():
		if minimumDistDict[ikey] <= miLevDist:
			subsetLevDistList.append(ikey)
		minimumDistArray.append(ikey + "\t"+ str(minimumDistDict[ikey]) +"\n")
	# for ifastaID in iFastaSeqDict.keys():
	# 	ifastaSeq = iFastaSeqDict[ifastaID]
	# 	mindist = 999
	# 	for jfastaID in iFastaSeqDict.keys():
	# 		jfastaSeq = iFastaSeqDict[jfastaID]
	# 		idist = -1
	# 		if not ifastaID == jfastaID:
	# 			idist = editdistance.eval(ifastaSeq, jfastaSeq)
	# 			if idist < mindist:
	# 				mindist = idist
	# 			if idist <= miLevDist: 
	# 				distArray.append(ifastaID +"\t"+ jfastaID +"\t"+ str(idist) +"\n")
	# 				subsetLevDistList.append(ifastaID)
	# 				subsetLevDistList.append(jfastaID)
	# 	if mindist <= miLevDist:
	# 		minimumDistArray.append(ifastaID + "\t"+ str(mindist) +"\n")
	# 	name2NameAccessionDict[ifastaID.split("/")[0]] = ifastaID
	return [name2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList]

## public function
def SeqSimilarityNetwork( miWorkDir, miQueryVersion, miOutputFolder, miSpecies, miLevDist):
	############################################
	## on a specific version, build Levenshtein distance matrix for each given species 
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miLevDist = int(miLevDist) ## if case user input a string

	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	versionFolder  = makeRightPath(miWorkDir+miQueryVersion) 

	## build sequence distance matrix
	for ispecies in miSpecies:
		matureLevDistList = []
		hairpinLevDistList = []
		childParentDictArray = []
		name2NameAccessionDict = {}
		for imolecular in ["mature", "hairpin"]:
			iFastaSeqDict = ReadFasta2DictSubsetKeepUFullKey((versionFolder+imolecular+".fa"), ispecies)
			[iname2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList] = calculatePairwiseLevenshteinDistance(iFastaSeqDict, imolecular,miLevDist)
			name2NameAccessionDict.update(iname2NameAccessionDict)
			## write out into file
			filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"+imolecular+"_"
			WriteList2File(distArray, filenamePrefix+"LevenshteinDistanceTable.tsv",miOutputFolder)
			WriteList2File(minimumDistArray, filenamePrefix+"minLevenshteinDistance.tsv",miOutputFolder)
			if imolecular == "mature":
				matureLevDistList = subsetLevDistList
			else:
				hairpinLevDistList = subsetLevDistList
		## build parent child relationship
		parentChildDict = ReadParentChildFromStr((versionFolder+"miRNA.str"), ispecies)

		matureLevDistUniqList = list(set(matureLevDistList))
		hairpinLevDistUniqList = list(set(hairpinLevDistList))
		for imature in matureLevDistUniqList:
			imatureID = imature.split("/")[0]
			imatureParents = parentChildDict[imatureID]
			# print imature +"\t"+ imatureParents
			for iparent in imatureParents.split(";"):
				iparentIDAcc = name2NameAccessionDict[iparent]
				if iparentIDAcc in hairpinLevDistUniqList:
					childParentDictArray.append(imature +"\t"+ iparentIDAcc +"\n")
		## write out into file
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		WriteList2File(childParentDictArray, filenamePrefix+"Mature2Hairpin.tsv",miOutputFolder)

		## ## ploting by R script --- MFEdistribution
		miRscriptsFolder = makeRightPath(miOutputFolder+"Rscripts")
		if not os.path.exists(miRscriptsFolder):
			os.makedirs(miRscriptsFolder)
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		assemblerRscript_SeqSimilarityNetwork(ispecies, filenamePrefix, miOutputFolder, miRscriptsFolder)
		rscriptExecuteCmd = " ".join(["Rscript" , (miRscriptsFolder+ispecies+"-SeqSimilarityNetwork.R")  ])
		os.system(rscriptExecuteCmd)

## public function
def BuildingCuratedGff( miWorkDir, miQueryVersion, miOutputFolder, miSpecies, miVersionList, miGffOptions):
	############################################
	## on a specific version, build curated annotation files for each given species 
	############################################
	parentalWorkDir = miWorkDir
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")
	
	## specify subfolder for each run of building a curated set
	miOutputFolder = makeRightPath(miOutputFolder+"curatedSet")
	if os.path.isdir(miOutputFolder):
		miOutputFolder = miOutputFolder.replace("curatedSet", "curatedSet.1")
		while os.path.isdir(miOutputFolder):
			index = miOutputFolder.split(".")[-1][:-1]
			miOutputFolder = makeRightPath(makeRightPath((makeRightPath(parentalWorkDir)+"resultsTables"))+"curatedSet."+str(int(index) +1))
		os.mkdir(miOutputFolder)
	else:
		os.mkdir(miOutputFolder)

	miSpecies = list(miSpecies)
	queryVersionFolder =  makeRightPath(miWorkDir+miQueryVersion) 
	for ispecies in miSpecies: ## iterate species
		excludingDict = {} ## excluding miRNA/pre-miRNA, and reasons
		
		## read in fasta into dictionary
		iMatureFastaDict = ReadFasta2DictSubsetKeepUFullKey((queryVersionFolder+"mature"+".fa"), ispecies)
		iHairpinFastaDict = ReadFasta2DictSubsetKeepUFullKey((queryVersionFolder+"hairpin"+".fa"), ispecies)
		
		## get keys from dictionary
		matureQueryList = iMatureFastaDict.keys()
		precursorQueryList = iHairpinFastaDict.keys()
		
		## child-parent convert, and name converting
		mirParentChildDict = ReadParentChildFromStr(queryVersionFolder+"miRNA.str", miSpecies)
		name2nameAccessionDict = buildName2NameAccessiondict(iHairpinFastaDict, iMatureFastaDict )
		
		## get profile from  overlapping	
		matureOverlappingArray = findingSequenceOverlapping(iMatureFastaDict, miQueryVersion, ispecies, "mature", mirParentChildDict, name2nameAccessionDict, iHairpinFastaDict)
		hairpinOverlappingArray = findingSequenceOverlapping(iHairpinFastaDict, miQueryVersion, ispecies, "hairpin", mirParentChildDict, name2nameAccessionDict, iHairpinFastaDict)
		
		## get profile from identical 
		uniqueMatureSeqDict = {}
		identicalSeqMatureArray = []
		uniqueHairpinSeqDict = {}
		identicalSeqHairpinArray = []
		for ifastaID in iMatureFastaDict.keys():
			ifastaSeq = iMatureFastaDict[ifastaID]
			if ifastaSeq in uniqueMatureSeqDict:
				uniqueMatureSeqDict[ifastaSeq] = uniqueMatureSeqDict[ifastaSeq] +";"+ ifastaID
			else:
				uniqueMatureSeqDict[ifastaSeq] = ifastaID
		for jfastaID in iHairpinFastaDict.keys():
			jfastaSeq = iHairpinFastaDict[jfastaID]
			if jfastaSeq in uniqueHairpinSeqDict:
				uniqueHairpinSeqDict[jfastaSeq] = uniqueHairpinSeqDict[jfastaSeq] +";"+ jfastaID
			else:
				uniqueHairpinSeqDict[jfastaSeq] = jfastaID
		for iseqkey in uniqueMatureSeqDict.keys():
			iseqvalue = uniqueMatureSeqDict[iseqkey]
			if ";" in iseqvalue:
				for i in iseqvalue.split(";"):
					identicalSeqMatureArray.append(i)
		for jseqkey in uniqueHairpinSeqDict.keys():
			jseqvalue = uniqueHairpinSeqDict[jseqkey]
			if ";" in jseqvalue:
				for j in jseqvalue.split(";"):
					identicalSeqHairpinArray.append(j)
		
		## get profile from "miRNA.diff"
		[ispeciesNEWlist, ispeciesNAMElist, ispeciesSEQUENCElist, ispeciesDELETElist] = readIDAClistFromDiffSubset((queryVersionFolder+"miRNA.diff"), ispecies)
		
		## get profile from Reverse Complementary
		[fullLenMatureRClist ,overlappingMatureRClist,  iReverseComplementaryMatureMiRArray] = findingReverseComplementary( iMatureFastaDict, miQueryVersion,ispecies, "mature")
		[fullLenHairpinRClist ,overlappingHairpinRClist,  iReverseComplementaryHairpinMiRArray] = findingReverseComplementary( iHairpinFastaDict, miQueryVersion,ispecies, "hairpin")

		querySpeciesGffPath = queryVersionFolder+"genomes/"+ispecies+".gff3"
		
		if not os.path.isfile(querySpeciesGffPath):
			querySpeciesGffPath = querySpeciesGffPath.replace("gff3", "gff2")
			if not os.path.isfile(querySpeciesGffPath):
				querySpeciesGffPath = querySpeciesGffPath.replace("gff2", "gff")
				if not os.path.isfile(querySpeciesGffPath):
					sys.exit("NO gff file found from your specification:"+querySpeciesGffPath)

		[headerLineArray, gffContentDict] = ReadGffContentKeyingID(querySpeciesGffPath,ispecies)
	
		print len(gffContentDict)
		## 
		miGffOptions = list(miGffOptions)
		for ioption in miGffOptions:
			## full length Reverse Complementary in miRNA
			if ioption == "RCm":
				unqualifyCount = 0
				for iquery in matureQueryList:
					if iquery in fullLenMatureRClist:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				imessage =  "in version %s,  %i pre-miRNAs are removed dut to full length reverse complementary in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
				print imessage
			
			## full length Reverse Complementary in pre-miRNA
			elif ioption == "RCp":
				for iquery in precursorQueryList:
					if iquery in fullLenHairpinRClist:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				imessage =  "in version %s,  %i pre-miRNAs are removed dut to full length reverse complementary in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
				print imessage
			
			## overlapping Reverse Complementary in miRNA
			elif ioption == "oRCm":
				for iquery in matureQueryList:
					if iquery in overlappingMatureRClist:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				imessage = "in version %s,  %i miRNAs are removed dut to overlapping reverse complementary in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
				print imessage
			
			## overlapping Reverse Complementary in pre-miRNA
			elif ioption == "oRCp":
				print "Reverse complementary in precursor sequence"
				for iquery in precursorQueryList:
					if iquery in overlappingHairpinRClist:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				imessage = "in version %s,  %i pre-miRNAs are removed dut to overlapping reverse complementary in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
				print imessage

			## pairwise Levenshtein Distance
			elif "LD" in ioption:
				maxdist = int(ioption.split("LD")[0])
				targetMolecular = ioption.split("LD")[1]
				if targetMolecular == "m":
					targetMolecular = "mature"
					[iname2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList] = calculatePairwiseLevenshteinDistance(iMatureFastaDict, targetMolecular, maxdist)
					unqualifyCount = 0
					for iquery in matureQueryList:
						if iquery in subsetLevDistList:
							unqualifyCount += 1
							excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					imessage = "in version %s,  %i miRNAs are removed dut to pairwise Levenshtein Distance less than %i in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, maxdist, ispecies)
					print imessage
				elif targetMolecular == "p":
					targetMolecular = "hairpin"
					[iname2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList] = calculatePairwiseLevenshteinDistance(iHairpinFastaDict, targetMolecular, maxdist)
					unqualifyCount = 0
					for iquery in precursorQueryList:
						if iquery in subsetLevDistList:
							unqualifyCount += 1
							excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					imessage = "in version %s,  %i pre-miRNAs are removed dut to pairwise Levenshtein Distance less than %i in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, maxdist, ispecies)
					print imessage
				elif targetMolecular == "mp" or targetMolecular == "pm":
					targetMolecular = ["mature", "hairpin"]
					[iname2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList] = calculatePairwiseLevenshteinDistance(iHairpinFastaDict, "hairpin", maxdist)
					unqualifyCount = 0
					for iquery in precursorQueryList:
						if iquery in subsetLevDistList:
							unqualifyCount += 1
							excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					imessage = "in version %s,  %i pre-miRNAs are removed dut to pairwise Levenshtein Distance less than %i in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, maxdist, ispecies)
					print imessage
					[iname2NameAccessionDict,  distArray,  minimumDistArray,  subsetLevDistList] = calculatePairwiseLevenshteinDistance(iMatureFastaDict, "mature", maxdist)
					unqualifyCount = 0
					for iquery in matureQueryList:
						if iquery in subsetLevDistList:
							unqualifyCount += 1
							excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					imessage = "in version %s,  %i miRNAs are removed dut to pairwise Levenshtein Distance less than %i in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, maxdist, ispecies)
					print imessage
				imessage = "setting maximum allowed pairwise Levenshtein Distance to %s for %s"%(maxdist, ";".join(targetMolecular))
				print imessage
			
			## version diff				
			elif ioption == "NAME":
				tagMatureSumCount = 0
				tagHairpinSumCount = 0
				for iquery in matureQueryList:
					if iquery in ispeciesNAMElist:
						tagMatureSumCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				for jquery in ispeciesNAMElist:
					if jquery in precursorQueryList:
						tagHairpinSumCount += 1
						excludingDict = UpdatingDict(jquery, ioption,excludingDict)
				print "in version %s, %i miRNAs and %i pre-miRNAs are removed dut to NAME change for the species %s"%(miQueryVersion, tagMatureSumCount, tagHairpinSumCount, ispecies)
			elif ioption == "NEW":
				tagMatureSumCount = 0
				tagHairpinSumCount = 0
				for iquery in matureQueryList:
					if iquery in ispeciesNEWlist:
						tagMatureSumCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				for jquery in ispeciesNEWlist:
					if jquery in precursorQueryList:
						tagHairpinSumCount += 1
						excludingDict = UpdatingDict(jquery, ioption,excludingDict)
				print "in version %s, %i miRNAs and %i pre-miRNAs are removed dut to NEW change for the species %s"%(miQueryVersion, tagMatureSumCount, tagHairpinSumCount, ispecies)
			elif ioption == "SEQUENCE":
				tagMatureSumCount = 0
				tagHairpinSumCount = 0
				for iquery in matureQueryList:
					if iquery in ispeciesSEQUENCElist:
						tagMatureSumCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				for jquery in ispeciesSEQUENCElist:
					if jquery in precursorQueryList:
						tagHairpinSumCount += 1
						excludingDict = UpdatingDict(jquery, ioption,excludingDict)
				print "in version %s, %i miRNAs and %i pre-miRNAs are removed dut to SEQUENCE change for the species %s"%(miQueryVersion, tagMatureSumCount, tagHairpinSumCount, ispecies)
			## ideentical and overlapping
			elif ioption == "ISm":
				unqualifyCount = 0
				for iquery in matureQueryList:
					if iquery in identicalSeqMatureArray:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "in version %s, %i pre-miRNAs  are removed dut to identical in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
			elif ioption == "ISp":
				unqualifyCount = 0
				for iquery in precursorQueryList:
					if iquery in identicalSeqHairpinArray:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "in version %s, %i pre-miRNAs  are removed dut to identical in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
			elif ioption == "OSm":
				unqualifyCount = 0
				for iquery in matureQueryList:
					if iquery in matureOverlappingArray:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "in version %s, %i miRNAs  are removed dut to overlapping in mature sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
			elif ioption == "OSp":
				unqualifyCount = 0
				for iquery in precursorQueryList:
					if iquery in hairpinOverlappingArray:
						unqualifyCount += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "in version %s, %i pre-miRNAs  are removed dut to overlapping in precursor sequence for the species %s"%(miQueryVersion, unqualifyCount, ispecies)
			## MFE
			elif "MFE" in ioption:
				minMFE = ioption.split("MFE")[0]
				maxMFE = ioption.split("MFE")[1]
				unqualifyMFEsum = 0
				iMFEhairpinDict = ReadMFEfromMiRStr(queryVersionFolder+"miRNA.str")
				for iquery in precursorQueryList:
					iqueryID = iquery.split("/")[0]
					if iqueryID in iMFEhairpinDict:
						iqueryMFE = iMFEhairpinDict[iqueryID]
						if iqueryMFE < minMFE or iqueryMFE > maxMFE:
							# precursorQueryList.remove(iquery)
							unqualifyMFEsum += 1
							excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					else:
						print "Could NOT find MFE for this precursor %s "%(iquery)
				print "in this species %s, %i precursors are un-qualify due to MFE limitation"%(ispecies, unqualifyMFEsum)
			## length range of mature sequence
			elif "ML" in ioption:
				minMatureLen = ioption.split("ML")[0]
				maxMatureLen = ioption.split("ML")[1]
				iMatureFastaDict = ReadFasta2DictSubsetKeepUFullKey((queryVersionFolder+"mature.fa"), ispecies)
				sumLenQualify = 0
				for iquery in  matureQueryList:
					if len(iMatureFastaDict[iquery]) < int(minMatureLen):
						sumLenQualify += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					elif  len(iMatureFastaDict[iquery])  > int(maxMatureLen):
						sumLenQualify += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "In species: %s,  %i miRNA  are excluded based on sequence length range %s to %s  in version %s"%(ispecies, sumLenQualify, minMatureLen, maxMatureLen, miQueryVersion)
			## length range of precursor sequence
			elif "PL" in ioption:
				minHairpinLen = ioption.split("PL")[0]
				maxHairpinLen = ioption.split("PL")[1]
				sumLenQualify = 0
				for iquery in  precursorQueryList:
					if len(iHairpinFastaDict[iquery]) < int(minHairpinLen):
						sumLenQualify += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
					elif  len(iHairpinFastaDict[iquery])  > int(maxHairpinLen):
						sumLenQualify += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				print "In species: %s,  %i miRNA precursors are excluded based on sequence length range %s to %s  in version %s"%(ispecies, sumLenQualify,minHairpinLen, maxHairpinLen, miQueryVersion)
			## genome coordinate
			elif ioption == "GCOR":
				locationDict = ReadCoordinate2Dict(queryVersionFolder+ "genomes/"+ispecies+".gff3")
				print len(locationDict)
				matureSumNoCoordinate = 0
				hairpinSumNoCoordinate = 0
				for iquery in matureQueryList:
					if not iquery in locationDict.keys():
						matureSumNoCoordinate += 1
						excludingDict = UpdatingDict(iquery, ioption,excludingDict)
				for jquery in precursorQueryList:
					if not jquery in locationDict.keys():
						hairpinSumNoCoordinate += 1
						excludingDict = UpdatingDict(jquery, ioption,excludingDict)
				print "In species: %s, %i miRNAs and %i miRNA precursors are missing genome coordinate in version %s"%(ispecies, matureSumNoCoordinate, hairpinSumNoCoordinate, miQueryVersion)
			## high confidence
			elif ioption == "HF":
				print "befor %i miRNAs and %i pre-miRNAs"%(len(matureQueryList), len(precursorQueryList))
				if not miQueryVersion in ["20", "21", "22"]:
					sys.exit("In your chosen miRBase version %s, no high confidence annotation available! Please choose another version."%(miQueryVersion))
				[highConfNameAccessionArray, highConfNameAccessionDict] = ReadAC2IDfromDatSubset(queryVersionFolder+ "miRNA_high_conf.dat", ispecies)
				print "high confidence %i entries in this species: %s"%(len(highConfNameAccessionArray), ispecies)
				nonHfMatureSum = 0
				for imature in matureQueryList:
					if not imature in highConfNameAccessionArray:
						nonHfMatureSum += 1
						excludingDict = UpdatingDict(imature, ioption,excludingDict)
				nonHfPrecursorSum = 0
				for ihairpin in precursorQueryList:
					if not ihairpin in highConfNameAccessionArray:
						nonHfPrecursorSum += 1
						excludingDict = UpdatingDict(ihairpin, ioption,excludingDict)
				print "In species: %s, %i miRNAs and %i miRNA precursors are removed from the list due to not labeled as high confidence in version %s"%(ispecies, nonHfMatureSum, nonHfPrecursorSum, miQueryVersion)
				print "after %i miRNAs and %i pre-miRNAs"%(len(matureQueryList), len(precursorQueryList))
			else:
				print "Unrecognised option: %s"%(ioption)
		## check the pre-miRNA of child mature miRNA still in the list of pre-miRNA
		## if not, remove this mature miRNA
		for imatureID in matureQueryList:
			if imatureID in excludingDict.keys():
				## excluding based on previous issues
				matureQueryList.remove(imatureID)
			else:
				## going to check if its parent pre-miRNA is good or not
				imatureIDParent = mirParentChildDict[imatureID.split("/")[0]]
				for iparent in imatureIDParent.split(";"):
					iparentAcc = name2nameAccessionDict[iparent]
					if not iparentAcc in precursorQueryList:
						if imatureID in matureQueryList:
							matureQueryList.remove(imatureID)
							excludingDict = UpdatingDict(imatureID, "unqualifyPrecursor",excludingDict)
		## check the child miRNA of pre-miRNA still in the list of miRNA
		## if not, remove this pre-miRNA
		for ihairpinID in precursorQueryList:
			if ihairpinID in excludingDict.keys():
				## excluding based on previous issues
				precursorQueryList.remove(ihairpinID)
			else:
				## going to check if its child miRNA is good or not
				ihairpinIDchild = mirParentChildDict[ihairpinID.split("/")[0]]
				for ichild in ihairpinIDchild.split(";"):
					ichildAcc = name2nameAccessionDict[ichild]
					if not ichildAcc in matureQueryList:
						if ihairpinID in precursorQueryList:
							precursorQueryList.remove(ihairpinID)
							excludingDict = UpdatingDict(ihairpinID, "unqualifyProduction",excludingDict)
		
		## build curated gff 
		gffOutputArray = []
		gffOutputArray += headerLineArray
		gffOutputArray.append("#Curated annotation criteria:\n#\t"+"; ".join(miGffOptions)+"\n#\n#\n")
		for ikey in matureQueryList:
			## in case some entries are missing coordinates, and GCOR is not included 
			if ikey in  gffContentDict.keys():
				gffOutputArray.append(gffContentDict[ikey])
			else:
				excludingDict = UpdatingDict(ikey, "GCOR",excludingDict)
				print "%s is annotated without genome coordinate"%(ikey)
		for jkey in precursorQueryList:
			if jkey in gffContentDict.keys():
				gffOutputArray.append(gffContentDict[jkey])
			else:
				excludingDict = UpdatingDict(jkey, "GCOR",excludingDict)
				print "%s is annotated without genome coordinate"%(jkey)
		gffOutputArray.append("\n")
		## write out into file
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		WriteList2File(gffOutputArray, filenamePrefix+"genome.gff3",miOutputFolder)

		## build reference mature fasta
		matureSeqOutputArray = []
		for imatureID in matureQueryList:
			imatureSeq = iMatureFastaDict[imatureID]
			matureSeqOutputArray.append(">"+imatureID+"\n"+imatureSeq.replace("U", "T")+"\n")
		## write out into file
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		WriteList2File(matureSeqOutputArray, filenamePrefix+"mature.fa",miOutputFolder)

		## build reference hairpin sequence 
		hairpinSeqOutputArray = []
		for ihairpinID in precursorQueryList:
			ihairpinSeq = iHairpinFastaDict[ihairpinID]
			hairpinSeqOutputArray.append(">"+ihairpinID+"\n"+ihairpinSeq.replace("U", "T")+"\n")
		## write out into file
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		WriteList2File(hairpinSeqOutputArray, filenamePrefix+"hairpin.fa",miOutputFolder)

		## output excluding, miRNA/pre-miRNA and their reasons
		excludingArray = []
		excludingArray.append("molecularID\treasons\n")
		for ikey in excludingDict.keys():
			excludingArray.append(ikey +"\t"+ excludingDict[ikey] + "\n")
		## write out into file
		filenamePrefix = "miRBase_V"+miQueryVersion +"_"+ispecies+"_"
		WriteList2File(excludingArray, filenamePrefix+"Excluding-miRNA-premiRNA.tsv",miOutputFolder)

def speciesValidation(miSpecies, miVersionList, miWorkDir):
	validatingMes = True
	for iversion in miVersionList:
		iVersionOrganismTxt = ReadOrganismToArray(makeRightPath(miWorkDir)+"miRBase/" +iversion +"/organisms.txt") 
		for ispecies in miSpecies:
			if ispecies not in iVersionOrganismTxt:
				if validatingMes:
					validatingMes = "%s is not "
				else:
					validatingMes += "; "


##  TMP demo test
def DemoTest( miWorkDir, miVersionList, miOutputFolder, miSpecies):
	############################################
	## tmp function
	############################################
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	miWorkDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	miRStrProfileArray = []

	for iversion in miVersionList: ## iterate versions
		versionStr = "V"+iversion
		iverionFolder = makeRightPath(miWorkDir+iversion) 

		for ispecies in miSpecies: ## iterate species
			iMirStr= ReadChild2ParentFromStr(iverionFolder+"miRNA.str", ispecies)
			multiplePrecursorsCount = 0
			for ikey in iMirStr.keys():
				if iMirStr[ikey].count(";") > 0:
					multiplePrecursorsCount += 1
				miRStrProfileArray.append(versionStr +"\t"+ ispecies +"\t"+ ikey +"\t"+ str(iMirStr[ikey].count(";")+1) +"\t"+ iMirStr[ikey] +"\n")
			print "in this species %s, %i have more than precursors out of %i"%(ispecies,multiplePrecursorsCount, len(iMirStr) )
	## write out into file
	filenamePrefix = "miRBase_"+"-".join(miVersionList) +"_"+"-".join(miSpecies)+"_"
	WriteList2File(miRStrProfileArray, filenamePrefix+"HairpinMoreThanMature-Profile.tsv",miOutputFolder)
	## 
	[highConfNameAccessionV22Array, highConfNameAccessionV22Dict] = ReadAC2IDfromDat(( makeRightPath(miWorkDir+"22")+ "miRNA_high_conf.dat"))
	queryFile = ["miRBase_V22_mmu_hairpin_minLevenshteinDistance.tsv", "miRBase_V22_hsa_hairpin_minLevenshteinDistance.tsv"]

	for ifile in queryFile:
		filehandle = open(miOutputFolder+ifile, "r")
		sumcount = 0
		print ifile
		for line in filehandle:
			iquery = line.split("\t")[0]
			if iquery in highConfNameAccessionV22Array:
				print iquery
				sumcount += 1
		filehandle.close()
		print "in file %s: %s moleculars are found in High Confidence set"%(ifile, str(sumcount))

## END for main program ##
####################################################


####################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## I/O function
def WriteList2File( yourList, yourFilename, yourFolder):
	############################################
	## write a list of string, output into file
	############################################
	putputfilename=""
	if yourFolder.endswith("/"):
		putputfilename = yourFolder + yourFilename
	else:
		putputfilename = yourFolder +"/"+ yourFilename
	filehandle = open(putputfilename, "w")
	filehandle.writelines(yourList)
	filehandle.close()

def ReadMFEfromMiRStr( queryMiRStrFile):
	############################################
	## read MFE information from miRNA.str file
	## return a dictionary, hairpin ID --- MFE
	############################################
	hairpin2MFEdict = {} ## hairpin name in key, MFE in value
	mirStrfile = open(queryMiRStrFile, "r")
	for line in mirStrfile:
		if line.startswith(">"):
			ihairpin = line.split()[0][1:]
			ispecies = ihairpin.split("-")[0]
			iMFEstr = line.split()[1][1:-1]
			hairpin2MFEdict[ihairpin] = str(iMFEstr)
	mirStrfile.close()
	return hairpin2MFEdict

def ReadParentChildFromStr( queryMiRStrFile, miSpecies):
	############################################
	## from a given miRNA.str file, read parental information between miRNA and pre-miRNA
	## return a dictionary: (1) miRNA ID --- pre-miRNA ID; (2) pre-miRNA ID --- miRNA ID
	############################################
	parentChildDict = {} ## 
	mirStrfile = open(queryMiRStrFile, "r")
	for line in mirStrfile:
		if line.startswith(">"):
			childList = []
			ihairpin = line.split()[0][1:]
			ispecies = ihairpin.split("-")[0]
			if ispecies in miSpecies:
				iMFEstr = line.split()[1][1:-1]
				for istring in line.split(":"):
					if "[" in istring:
						childList.append(istring.split("[")[-1]) 
				parentChildDict[ihairpin] = ";".join(childList)
				for ichild in childList:
					if ichild in parentChildDict.keys():
						parentChildDict[ichild] = parentChildDict[ichild]+";"+ ihairpin
					else:
						parentChildDict[ichild] =  ihairpin
	mirStrfile.close()
	return parentChildDict

def ReadChild2ParentFromStr( queryMiRStrFile, miSpecies):
	############################################
	## from a given miRNA.str file, read parental information between miRNA and pre-miRNA
	## return a dictionary, ONLY for: (1) miRNA ID --- pre-miRNA ID 
	############################################
	parentChildDict = {} ## 
	mirStrfile = open(queryMiRStrFile, "r")
	for line in mirStrfile:
		if line.startswith(">"):
			childList = []
			ihairpin = line.split()[0][1:]
			ispecies = ihairpin.split("-")[0]
			if ispecies in miSpecies:
				iMFEstr = line.split()[1][1:-1]
				for istring in line.split(":"):
					if "[" in istring:
						childList.append(istring.split("[")[-1]) 
				for ichild in childList:
					if ichild in parentChildDict.keys():
						parentChildDict[ichild] = parentChildDict[ichild]+";"+ ihairpin
					else:
						parentChildDict[ichild] =  ihairpin
	mirStrfile.close()
	return parentChildDict

def ReadFasta2Dict( fastaFilename):
	############################################
	## Read in a fasta file, return a dictionary, key --  name, value -- sequence .  also subset by a given species. replace the U with T
	############################################
	fastaDict = {}
	fastaHandle = open(fastaFilename, "r")
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		fastaDict[ifasta.id] = str(ifasta.seq).replace("U", "T")
	fastaHandle.close()
	print "read in %i sequences from %s "%(len(fastaDict), fastaFilename)
	return fastaDict

def ReadFasta2DictFeaturesSubset( fastaFilename, speciesCode):
	############################################
	## Read in a fasta file, return a dictionary, key --  name, value -- features including accession number, sequence length, GC content, 
	############################################
	fastaDict = {}
	fastaHandle = open(fastaFilename, "r")
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		if (ifasta.id).split("-")[0] == speciesCode:
			iAccession = (ifasta.description).split()[1]
			iSeqLength = len(str(ifasta.seq))
			igccontent = GetGCcontent(str(ifasta.seq))
			fastaDict[ifasta.id] = str(iAccession) +"\t"+ str(iSeqLength) +"\t"+ str(igccontent) 
	fastaHandle.close()
	print "read in %i sequences from %s for %s"%(len(fastaDict), fastaFilename, speciesCode)
	return fastaDict

def ReadFasta2DictSubsetKeepUFullKey( fastaFilename, speciesCode):
	############################################
	## Read in a fasta file, return a dictionary, key --  name/accession ID, value -- sequence .  also subset by a given species. keep the U
	############################################
	fastaDict = {}
	fastaHandle = open(fastaFilename, "r")
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		if (ifasta.id).split("-")[0] ==  speciesCode  :
			iaccession = (ifasta.description).split()[1]
			fastaDict[ifasta.id+"/"+iaccession] = str(ifasta.seq)
	fastaHandle.close()
	print "read in %i sequences from %s in the %s species"%(len(fastaDict), fastaFilename, speciesCode)
	return fastaDict

def ReadFasta2DictSubsetFullKey( fastaFilename, speciesCode):
	############################################
	## Read in a fasta file, return a dictionary, key --  name/accession ID, value -- sequence .  also subset by a given species. keep the U
	############################################
	fastaDict = {}
	fastaHandle = open(fastaFilename, "r")
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		if (ifasta.id).split("-")[0] ==  speciesCode  :
			iaccession = (ifasta.description).split()[1]
			fastaDict[ifasta.id+"/"+iaccession] = str(ifasta.seq).replace("U", "T")
	fastaHandle.close()
	print "read in %i sequences from %s in the %s species"%(len(fastaDict), fastaFilename, speciesCode)
	return fastaDict

def ReadFasta2DictFullKey( fastaFilename):
	############################################
	## Read in a fasta file, return a dictionary, key --  name/accession ID, value -- sequence .  also subset by a given species. keep the U
	############################################
	fastaDict = {}
	fastaHandle = open(fastaFilename, "r")
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		# if (ifasta.id).split("-")[0] ==  speciesCode  :
		iaccession = (ifasta.description).split()[1]
		fastaDict[ifasta.id+"/"+iaccession] = str(ifasta.seq).replace("U", "T")
	fastaHandle.close()
	# print "read in %i sequences from %s in the %s species"%(len(fastaDict), fastaFilename, speciesCode)
	return fastaDict

def ReadAC2IDfromDat( datFileName):
	############################################
	## read the information about accession-to-ID/Accession
	############################################
	AC2IDdictionary = {}
	nameAcccessionArray = []
	datHandle  = open(datFileName,"r")
	entryID = "NA"
	entryAC = "NA"
	for line in datHandle:
		line = line.rstrip()
		if line.startswith("ID"):
			entryID = line[5:].split()[0]
		elif line.startswith("AC"):
			entryAC = line[5:-1]
			if entryAC.count(";") > 1:
				print entryAC
				for i in entryAC.split(";"):
					print entryID +"/"+ i
			else:
				nameAcccessionArray.append(entryID +"/"+ entryAC)
				AC2IDdictionary[entryAC] = entryID +"/"+ entryAC
	datHandle.close()
	return [nameAcccessionArray, AC2IDdictionary]

def ReadAC2IDfromDat_Full( datFileName):
	############################################
	## read the information about accession-to-ID/Accession
	## also including mature miRNA accession
	############################################
	precursorAC2IDdict = {}
	precursorAccArray = []
	matureAC2IDdict = {}
	matureAccArray = []
	datHandle  = open(datFileName,"r")
	precursorEntryID = "NA"
	precursorEntryAC = "NA"
	matureEntryID = "NA"
	matureEntryAC = "NA"
	for line in datHandle:
		line = line.rstrip()
		if line.startswith("ID"):
			precursorEntryID = line[5:].split()[0]
		elif line.startswith("AC"):
			precursorEntryAC = line[5:-1]
			if precursorEntryAC.count(";") > 1:
				print precursorEntryAC
				for i in precursorEntryAC.split(";"):
					print precursorEntryID +"/"+ i
			else:
				precursorAccArray.append(precursorEntryAC)
				precursorAC2IDdict[precursorEntryAC] = precursorEntryID +"/"+ precursorEntryAC

		elif line.startswith("FT"):
			if "/accession=" in line:
				matureEntryAC = line.split("=")[1].replace("\"", "")
			elif "/product=" in line:
				matureEntryID = line.split("=")[1].replace("\"", "")
				matureAccArray.append(matureEntryAC)
				matureAC2IDdict[matureEntryAC] = matureEntryID +"/"+ matureEntryAC
	datHandle.close()
	return [precursorAccArray, precursorAC2IDdict, matureAccArray, matureAC2IDdict]

def ReadAC2IDfromDatSubset( datFileName, speciesCode):
	############################################
	## read the information about accession-to-ID/Accession
	## subset for species
	############################################
	AC2IDdictionary = {}
	nameAcccessionArray = []
	parent2childMirDict = {}
	datHandle  = open(datFileName,"r")
	entryHairpinID = "NA" ##  parental 
	entryHairpinAC = "NA"
	entryMatureID = "NA"
	entryMatureAC = "NA"
	for line in datHandle:
		line = line.rstrip()
		if line.startswith("ID"):
			entryHairpinID = line[5:].split()[0]
		elif line.startswith("AC"):
			entryHairpinAC = line[5:-1]
			if entryHairpinAC.count(";") > 1:
				print entryHairpinAC
				for i in entryHairpinAC.split(";"):
					print entryHairpinID +"/"+ i
			else:
				if entryHairpinID.split("-")[0] == speciesCode:
					nameAcccessionArray.append(entryHairpinID +"/"+ entryHairpinAC)
					AC2IDdictionary[entryHairpinAC] = entryHairpinID +"/"+ entryHairpinAC
		elif line.startswith("FT"):
			if "/accession=" in line:
				entryMatureAC = line.split("=")[1][1:-1]
			elif "/product=" in line:
				entryMatureID = line.split("=")[1][1:-1]
				if entryMatureID.split("-")[0] == speciesCode:
					nameAcccessionArray.append(entryMatureID +"/"+ entryMatureAC)
					AC2IDdictionary[entryHairpinAC] = entryMatureID +"/"+ entryMatureAC
					if (entryHairpinID +"/"+ entryHairpinAC) in parent2childMirDict.keys():
						parent2childMirDict[entryHairpinID +"/"+ entryHairpinAC] = parent2childMirDict[entryHairpinID +"/"+ entryHairpinAC] + entryMatureID +"/"+ entryMatureAC
					else:
						parent2childMirDict[entryHairpinID +"/"+ entryHairpinAC] =  entryMatureID +"/"+ entryMatureAC	
	datHandle.close()
	return [nameAcccessionArray, AC2IDdictionary]

def ReadCoordinate2Dict( gffFilename):
	############################################
	## reading genome location of hairpin and mature both
	############################################
	locationDict = {}
	gffHandle = open(gffFilename,"r")
	for iline in gffHandle:
		if not iline.startswith("#"):
			iline = iline.rstrip()
			dtype = iline.split("\t")[2]
			if dtype == "miRNA_primary_transcript":
				chroStr = iline.split("\t")[0]
				posirionStr = str(iline.split("\t")[3]) +"-" +str(iline.split("\t")[4])
				strandStr = iline.split("\t")[6]
				jfeatures = iline.split("\t")[8].split(";")
				hairpinName = jfeatures[2].split("=")[1]
				hairpinID = jfeatures[1].split("=")[1]
				hairpinMergedName = hairpinName+"/"+hairpinID
				if hairpinMergedName in locationDict:
					locationDict[hairpinMergedName] = locationDict[hairpinMergedName]+ ";"+ chroStr +":"+posirionStr +":"+strandStr
				else:
					locationDict[hairpinMergedName] = chroStr +":"+posirionStr +":"+strandStr
			elif dtype == "miRNA":
				chroStr = iline.split("\t")[0]
				posirionStr = str(iline.split("\t")[3]) +"-" +str(iline.split("\t")[4])
				strandStr = iline.split("\t")[6]
				jfeatures = iline.split("\t")[8].split(";")
				hairpinName = jfeatures[2].split("=")[1]
				hairpinID = jfeatures[1].split("=")[1]
				hairpinMergedName = hairpinName+"/"+hairpinID
				if hairpinMergedName in locationDict:
					locationDict[hairpinMergedName] = locationDict[hairpinMergedName]+ ";"+ chroStr +":"+posirionStr +":"+strandStr
				else:
					locationDict[hairpinMergedName] = chroStr +":"+posirionStr +":"+strandStr
	gffHandle.close()
	return locationDict

def GetGCcontent( querySeqStr):
	############################################
	## given a sequence string, calculate the GC content of it; return GC content without percentage
	############################################
	querySeqStr = querySeqStr.upper()
	countingGC = 0
	countingGC += querySeqStr.count("G") 
	countingGC += querySeqStr.count("C")
	gccontent = float(countingGC)/len(querySeqStr)
	return str(gccontent)

def readIDAClistFromDiffSubset( queryDiffFilePath, speciesCode):
	############################################
	## read profile from "miRNA.diff" file
	############################################
	list4NEW = []
	list4NAME = []
	list4SEQUENCE = []
	list4DELETE = []

	diffHandler = open(queryDiffFilePath, "r")
	for line in diffHandler:
		if not line.startswith("#"):
			line = line.rstrip()
			lineSplitted = line.split("\t")
			ispecies = lineSplitted[1].split("-")[0]
			if ispecies == speciesCode:  ## select the species
				mergedIDname = lineSplitted[1] + "/"+ lineSplitted[0]
				for i in range(2,len(lineSplitted)): ## in case, some entry can has more than one tags, like NAME SEQUENCE
					itag = lineSplitted[i]
					if itag == "NEW":
						list4NEW.append(mergedIDname)
					elif itag == "NAME":
						list4NAME.append(mergedIDname)
					elif itag == "SEQUENCE":
						list4SEQUENCE.append(mergedIDname)
					else :
						list4DELETE.append(mergedIDname)
	diffHandler.close()
	return [list4NEW, list4NAME, list4SEQUENCE, list4DELETE]

def ReadGffV3GetGenomeCoordinates( queryGffPath):
	############################################
	## Read genome coordinate for miRNA/pre-miRNA from gff3
	############################################
	gffHandle = open(queryGffPath,"r")
	matureEntryCount = 0
	hairpinEntryCount = 0
	matureMultiLocationDict = {}
	hairpinMultiLocationDict = {}
	for line in gffHandle:
		if not line.startswith("#"): ## skip header
			line = line.rstrip()
			lineSplitted = line.split("\t")
			molecularType = lineSplitted[2]
			if molecularType == "miRNA": ## miRNA entry
				matureEntryCount += 1
				ifeatures =lineSplitted[8].split(";")
				matureName = ifeatures[2].split("=")[1]
				matureAccessionID = ifeatures[0].split("=")[1]
				matureAliasID = ifeatures[1].split("=")[1]
				hairpinID = ifeatures[3].split("=")[1]
				genemoCoordinateStr = lineSplitted[0] + ":"+lineSplitted[3] + "-"+lineSplitted[4]+ ":"+lineSplitted[6] + "@"+hairpinID
				fullMatureName = matureName+"/"+matureAliasID
				if fullMatureName in matureMultiLocationDict:
					matureMultiLocationDict[fullMatureName] = matureMultiLocationDict[fullMatureName] +";" + genemoCoordinateStr
				else:
					matureMultiLocationDict[fullMatureName] = genemoCoordinateStr
			elif molecularType == "miRNA_primary_transcript":	## pre-miRNA entry
				hairpinEntryCount += 1
				jfeatures = lineSplitted[8].split(";")
				hairpinName = jfeatures[2].split("=")[1]
				hairpinAccessionID = jfeatures[0].split("=")[1]
				hairpinAliasID = jfeatures[1].split("=")[1]
				fullHairpinName = hairpinName+"/"+hairpinAliasID 
				genemoCoordinateStr = lineSplitted[0] + ":"+lineSplitted[3] + "-"+lineSplitted[4]+ ":"+lineSplitted[6] + "@"+hairpinAliasID
				if fullHairpinName in hairpinMultiLocationDict:
					hairpinMultiLocationDict[fullHairpinName] = hairpinMultiLocationDict[fullHairpinName] +";"+genemoCoordinateStr
				else:
					hairpinMultiLocationDict[fullHairpinName] = genemoCoordinateStr
	gffHandle.close()
	return [ matureEntryCount,  hairpinEntryCount,  matureMultiLocationDict, hairpinMultiLocationDict]

def ReadGffV2GetGenomeCoordinates( queryGffPath):
	############################################
	## Read genome coordinate for miRNA/pre-miRNA from gff2 or gff
	## gff2/gff does NOT have information pre-miRNA
	############################################
	gffHandle = open(queryGffPath,"r")
	matureEntryCount = 0
	hairpinEntryCount = 0
	matureMultiLocationDict = {}
	hairpinMultiLocationDict = {}
	for line in gffHandle:
		if not line.startswith("#"): ## skip header
			line = line.rstrip()
			lineSplitted = line.split("\t")
			molecularType = lineSplitted[2]
			if molecularType == "miRNA": ## miRNA entry
				matureEntryCount += 1
				ifeatures =lineSplitted[8].split(";")
				matureName = ifeatures[1].split("=")[1][1:-1]
				matureAccessionID = ifeatures[0].split("=")[1][1:-1]
				genemoCoordinateStr = lineSplitted[0] + ":"+lineSplitted[3] + "-"+lineSplitted[4]+ ":"+lineSplitted[6] 
				fullMatureName = matureName+"/"+matureAccessionID
				if fullMatureName in matureMultiLocationDict:
					matureMultiLocationDict[fullMatureName] = matureMultiLocationDict[fullMatureName] +";" + genemoCoordinateStr
				else:
					matureMultiLocationDict[fullMatureName] = genemoCoordinateStr

			elif molecularType == "miRNA_primary_transcript":	## pre-miRNA entry, actually seems in gff2 and gff there is not annotated genome coordinates for pre-miRNA
				hairpinEntryCount += 1

				jfeatures = lineSplitted[8].split(";")
				hairpinName = jfeatures[1].split("=")[1][1:-1]
				hairpinAccessionID = jfeatures[0].split("=")[1][1:-1]
				# hairpinAliasID = jfeatures[1].split("=")[1][1:-1]
				fullHairpinName = hairpinName+"/"+hairpinAccessionID 
				genemoCoordinateStr = lineSplitted[0] + ":"+lineSplitted[3] + "-"+lineSplitted[4]+ ":"+lineSplitted[6] 
				#+ "@"+hairpinAliasID
				if fullHairpinName in hairpinMultiLocationDict:
					hairpinMultiLocationDict[fullHairpinName] = hairpinMultiLocationDict[fullHairpinName] +";"+genemoCoordinateStr
				else:
					hairpinMultiLocationDict[fullHairpinName] = genemoCoordinateStr

	gffHandle.close()
	return [ matureEntryCount,  hairpinEntryCount,  matureMultiLocationDict, hairpinMultiLocationDict]

def UpdatingDict( iqueryKey, iqueryValue, queryDict):
	############################################
	## update the queryDict with iqueryKey and iqueryValue
	############################################
	if iqueryKey in queryDict.keys():
		queryDict[iqueryKey] = queryDict[iqueryKey] + ";" +iqueryValue
	else:
		queryDict[iqueryKey] = iqueryValue
	return queryDict

def ReadGffContentKeyingID( querySpeciesGffPath, querySpecies):
	############################################
	## read gff content from a given gff file for a species
	## return a header line list and a dictionary of ID key with gff content value 
	## not care it's miRNA or pre-miRNA
	############################################
	gffHandler = open(querySpeciesGffPath, "r")
	##  start process the gff file 
	headerLineArray = []
	gffContentDict = {}
	if querySpeciesGffPath.endswith("gff3"):
		## gff3 format
		for line in gffHandler:
			if line.startswith("#"):
				headerLineArray.append(line)
			else:
				lineSplitted = (line.rstrip()).split("\t")
				molecularType = lineSplitted[2]
				ifeatures = lineSplitted[8].split(";")
				molecularID = ifeatures[2].split("=")[1]
				thisSpecies = molecularID.split("-")[0]
				molecularAlias = ifeatures[1].split("=")[1]
				if thisSpecies == querySpecies:
					mergedName = molecularID+"/"+molecularAlias
					if mergedName in gffContentDict.keys():
						gffContentDict[mergedName] = gffContentDict[mergedName] + line
					else:
						gffContentDict[mergedName] = line
	else:
		## gff2 or gff format
		for line in gffHandler:
			if line.startswith("#"):
				headerLineArray.append(line)
			else:
				lineSplitted = line.split("\t")
				molecularType = lineSplitted[2]
				ifeatures =lineSplitted[8].split(";")
				molecularID = ifeatures[1].split("=")[1][1:-1]
				thisSpecies = molecularID.split("-")[0]
				molecularAccession = ifeatures[0].split("=")[1][1:-1]
				if thisSpecies == querySpecies:
					mergedName = molecularID+"/"+molecularAccession
					if mergedName in gffContentDict.keys():
							gffContentDict[mergedName] = gffContentDict[mergedName] + line
					else:
						gffContentDict[mergedName] = line
	gffHandler.close()
	return [headerLineArray, gffContentDict]

def ReadChildParent( queryMirStrFile):
	############################################
	## read child and parent relationship from miRNA.str
	## and counting miRNAs and pre-miRNAs in each species
	############################################
	mirstrHandle = open(queryMirStrFile, "r")

	childParentSummary = {}

	AllChildlist = []
	AllParentlist = []

	for line in mirstrHandle:
		if line.startswith(">"):
			line = line.rstrip()
			lineSplitted = line.split(" ")
			hairpinID = lineSplitted[0][1:]
			AllParentlist.append( hairpinID )
			ispecies = hairpinID.split("-")[0]
			childList = []
			for i in range(4,len(lineSplitted)):
				childList.append( lineSplitted[i][1:-1].split(":")[0])
				AllChildlist.append(lineSplitted[i][1:-1].split(":")[0])

			singleProduction = 0
			if len(childList) > 1:
				singleProduction = 1
			if ispecies in childParentSummary.keys():
				tmplist = childParentSummary[ispecies]
				tmplist[0] += 1 ##  the number of pre-miRNAs
				tmplist[1] += len(childList) ## the number of miRNAs
				tmplist[2] += singleProduction ## the number of pre-miRNAs that have both 5p/3p
				childParentSummary[ispecies] = tmplist
			else:
				tmplist = [1,len(childList), singleProduction]
				childParentSummary[ispecies] = tmplist
	# childParentSummaryArrayDict =\
	return [childParentSummary, AllChildlist, AllParentlist]

def ReadFastaNameAsArray( yourFastaFile):
	fastaIdArray = []
	fastaHandle = open(yourFastaFile)
	for irecord in SeqIO.parse(fastaHandle,"fasta"):
		fastaIdArray.append(irecord.id)
	fastaHandle.close()
	return fastaIdArray

def BuildFastaLengthDict( yourFastaFile):
	fastaLenDict = {}
	fastaHandle = open(yourFastaFile)
	for ifasta in SeqIO.parse(fastaHandle,"fasta"):
		if not ifasta.id in fastaLenDict:
			fastaLenDict[ifasta.id] = len(ifasta.seq)
		else:
			print "should not have duplicate name"
			sys.exit()
	return fastaLenDict

def MFEdensityInMiRStr( miWorkDir, miVersionList, miOutputFolder):
	############################################
	## read MFE information from miRNA.str file
	## build a MFE density/histgram for each species
	############################################
	## prepare output folder and work directory
	if os.path.isdir(miOutputFolder):
		miOutputFolder = makeRightPath(miOutputFolder)
	else:
		miOutputFolder = makeRightPath(miWorkDir)+"resultsTables"
		miOutputFolder = makeRightPath(miOutputFolder)		
	if not os.path.exists(miOutputFolder):
		os.makedirs(miOutputFolder)
	mirbaseDir = makeRightPath(makeRightPath(miWorkDir) +"miRBase")

	speciesMFEarray = [] 

	for iversion in miVersionList:
		iMirStrFile = makeRightPath(mirbaseDir+iversion)+"miRNA.str"
		if not os.path.isfile(iMirStrFile):
			print "SKIPPING! Could not find miRNA.str file in version %s for processing: %s."%(iversion, iMirStrFile)

		speciesDecodingDict = ReadOrganismDecoding( makeRightPath(mirbaseDir+iversion)+"organisms.txt")
		mirStrfile = open(iMirStrFile, "r")
		for line in mirStrfile:
			if line.startswith(">"):
				ihairpin = line.split()[0][1:]
				ispecies = ihairpin.split("-")[0]
				iMFEstr = line.split()[1][1:-1]
				# print speciesDecodingDict[ispecies][1]
				speciesMFEarray.append("V"+iversion+"\t"+ ispecies +"\t"+ speciesDecodingDict[ispecies][1] +"\t"+ speciesDecodingDict[ispecies][2]+"\t"+ iMFEstr +"\n")
				
		mirStrfile.close()
	filenamePrefix = "miRBase_"+"V"+str(miVersionList[0]) + "-"+"V"+str(miVersionList[-1]) +"_"
	WriteList2File(speciesMFEarray, filenamePrefix+"MFE4SpeciesInVersions.tsv",miOutputFolder)

def ReadOrganismDecoding( queryOrganismTxt):
	############################################
	## read in orhanism.txt to decode the species code 
	############################################
	if not os.path.isfile(queryOrganismTxt):
		sys.exit("Could not find organism.txt file!")
	txtHandle = open(queryOrganismTxt)
	speciesDecodingDict = {}
	for line in txtHandle:
		if not line.startswith("#"):
			line = line.rstrip()
			lineSplitted = line.split("\t")
			organismFeatures = []
			organismFeatures.append(lineSplitted[2])
			speciesDecodingDict[lineSplitted[0]] = organismFeatures + lineSplitted[3].split(";")

	txtHandle.close()
	return speciesDecodingDict

def ReadOrganismToArray( queryOrganismTxt):
	############################################
	## read in orhanism.txt to decode the species code 
	############################################
	if not os.path.isfile(queryOrganismTxt):
		sys.exit("Could not find organism.txt file!")
	txtHandle = open(queryOrganismTxt)
	speciesArray = {}
	for line in txtHandle:
		if not line.startswith("#"):
			line = line.rstrip()
			lineSplitted = line.split("\t")
			organismFeatures = []
			speciesArray.append(lineSplitted[0] ) 
	txtHandle.close()
	return speciesArray
	

## END ##
####################################################


####################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## R script assembling
def assemblerRscript_countingMatureInVersion( filenamePrefix, miOutputFolder, miRscriptsFolder):
	# filenamePrefix, miOutputFolder, miRscriptsFolder
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"countingMatureInSpeciesInVersion.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("speciesCounting <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(speciesCounting) <- c('version', 'species',  'count')")
	scriptArray.append("## set version factor")
	scriptArray.append("speciesCounting$version <- factor(speciesCounting$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 20, height = 10)")
	scriptArray.append("ggplot(speciesCounting, aes(x=species , y=version,fill=count)) + geom_tile(stat = 'identity') +scale_fill_gradient2(low='green',  mid='yellow',high= 'red', midpoint=500, name='mature count', na.value = 'grey50', guide = 'colourbar') + theme(axis.text.x = element_text(size = 5,hjust = 0,vjust = 0,angle = 270), legend.position = 'bottom' ,axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.text=element_text(size=15,angle = 270),legend.title=element_text(size=18)) ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	inputfile = miOutputFolder+filenamePrefix+"countingSpeciesInVersion.tsv"
	scriptArray.append("#")
	scriptArray.append("## read in data")
	scriptArray.append("speciesCountingSum <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("colnames(speciesCountingSum) <- c(\"version\", \"count\")")
	scriptArray.append("speciesCountingSum$version <- factor(speciesCountingSum$version,  levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append("ggplot(speciesCountingSum,aes(y=count,x=version)) + geom_bar(stat=\"identity\")  +theme(axis.text.x = element_text(size = 12,hjust = 0,vjust = 0), axis.text.y = element_text(size = 12), axis.title.x = element_text(size=20), axis.title.y = element_text(size=15) ) +xlab(\"Versions\") +ylab(\"Species count\")")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	inputfile = miOutputFolder+filenamePrefix+"countingMatureInVersion.tsv"
	scriptArray.append("#")
	scriptArray.append("## read in data")
	scriptArray.append("matureCountingSum <- read.delim(\""+inputfile+"\", sep='\t',header=TRUE)")
	scriptArray.append("colnames(matureCountingSum) <- c(\"version\", \"count\")")
	scriptArray.append("matureCountingSum$version <- factor(matureCountingSum$version,  levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append("ggplot(matureCountingSum,aes(y=count,x=version)) + geom_bar(stat = \"identity\")  +theme(axis.text.x = element_text(size = 12,hjust = 0,vjust = 0), axis.text.y = element_text(size = 12) ,axis.title.x = element_text(size=20),axis.title.y = element_text(size=15) ) +xlab(\"Versions\") +ylab(\"Species count\") ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "countingMatureInVersion.R",  miRscriptsFolder)

def assemblerRscript_matureSeqLenDistribution( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"MatureSeqLenDistInVersions.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("matureLenDistTable <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(matureLenDistTable) <- c(\"version\",\"species\", \"len\", \"count\")")
	scriptArray.append("## set version factor")
	scriptArray.append("matureLenDistTable$version <- factor(matureLenDistTable$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("## sum rows by species")
	scriptArray.append("speciesTotal <- aggregate(count~species, data=matureLenDistTable, FUN=sum)")
	scriptArray.append("## order species by sum")
	scriptArray.append("speciesTotal <- speciesTotal[order(speciesTotal$count),] ")
	scriptArray.append("speciesGR500 <- speciesTotal[which(speciesTotal$count > 500),]$species")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 20, height = 10)")
	scriptArray.append(" ggplot(matureLenDistTable[which(matureLenDistTable$species %in% speciesGR500),], aes(y=len, x=species, size=count))  +geom_point() +ylab(\"Sequence length\") +xlab(\"Species\") +theme(axis.text.x = element_text(size = 9,hjust = 0,vjust = 0.5,angle = 270), axis.title.x = element_text(size=20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size=20) ) + facet_wrap(~version)  ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "matureSeqLenDistribution.R",  miRscriptsFolder)

def assemblerRscript_buildText4SpeciesCloud( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"Text4SpeciesCloud.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('tm')) install.packages('tm')")
	scriptArray.append("library('tm')")
	scriptArray.append("if (!require('wordcloud')) install.packages('wordcloud')")
	scriptArray.append("library('wordcloud')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("speciesAbundText <- readLines(\""+inputfile+"\" )")
	scriptArray.append("")
	scriptArray.append("speciesAbundCorpus = Corpus(VectorSource(speciesAbundText))")
	scriptArray.append("speciesAbundCorpus = tm_map(speciesAbundCorpus, tolower)")
	scriptArray.append("speciesAbundCorpus = tm_map(speciesAbundCorpus, removePunctuation)")
	scriptArray.append("speciesAbundCorpus = tm_map(speciesAbundCorpus, removeNumbers)")
	scriptArray.append("speciesAbundCorpus = tm_map(speciesAbundCorpus, removeWords, stopwords(\"english\"))")
	scriptArray.append("myDTM = TermDocumentMatrix(speciesAbundCorpus, control = list(minWordLength = 1))")
	scriptArray.append("m = as.matrix(myDTM)")
	scriptArray.append("v = sort(rowSums(m), decreasing = TRUE)")
	scriptArray.append("")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 4, height = 4)")
	scriptArray.append("set.seed(4363)  ")
	scriptArray.append("wordcloud(names(v), v, min.freq = 1) ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "buildText4SpeciesCloud.R",  miRscriptsFolder)

def assemblerRscript_diffCountMatureHairpin( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"CountingMatureHairpinBySpeciesInVersions.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("if (!require('scales')) install.packages('scales')")
	scriptArray.append("library('scales')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("miRBase.inVersions.countingMatureHairpin <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(miRBase.inVersions.countingMatureHairpin) <- c(\"version\",\"species\", \"matureCount\", \"hairpinCount\",\"diff\")")
	scriptArray.append("## set version factor")
	scriptArray.append("miRBase.inVersions.countingMatureHairpin$version <- factor(miRBase.inVersions.countingMatureHairpin$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 20, height = 10)")
	scriptArray.append(" ggplot(miRBase.inVersions.countingMatureHairpin, aes(y= matureCount, x= hairpinCount, color = diff, label= species)) +geom_point() + geom_abline(intercept = 0, slope = 1, color=\"darkgreen\", linetype=\"dashed\", size=0.5,alpha=0.7) +theme(axis.text.x = element_text(size=10,hjust = 1,vjust = 0),axis.text.y = element_text(size=10),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))  +scale_color_gradientn(colours=c(\"blue\",\"cyan\",\"green\", \"yellow\",\"red\") ,values=rescale(c(-80,-40,0,100,800)))   +facet_wrap(~version,nrow = 2) +ylab(\"mature count\")+xlab(\"hairpin count\")  ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "diffCountMatureHairpin.R",  miRscriptsFolder)

def assemblerRscript_MFEdistribution( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"fullDistribution4MFE.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("if (!require('scales')) install.packages('scales')")
	scriptArray.append("library('scales')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("mfeDistTable <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(mfeDistTable) <- c(\"version\",\"species\", \"precuror\", \"AccessionID\",\"MFE\",\"hairpinLength\",\"GCcontent\",\"adjMFE\",\"MFEIndex\" )")
	scriptArray.append("## set version factor")
	scriptArray.append("mfeDistTable$version <- factor(mfeDistTable$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append(" ggplot(mfeDistTable[which(mfeDistTable$species %in% c(\"hsa\",\"mmu\")),], aes(x=version, y=MFE))+ geom_violin(scale = \"count\", stat = \"ydensity\", color=\"darkblue\", width=1, fill=\"cyan\")+ geom_jitter(width = 0.4, alpha=0.2, size=0.5) + stat_summary(fun.y = mean, geom=\"point\",color=\"orange\")+ facet_wrap(~species,ncol=1)+ ylab(\"Estimated Minimum Free Energy of hairpin\") +xlab(\"Version\") + theme(strip.text.x = element_text(size = 25), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +geom_abline(intercept = -37, slope = 0, color=\"red\", linetype=\"dashed\", size=1) ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "MFEdistribution.R",  miRscriptsFolder)

def assemblerRscript_identicalWithinSpecies( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"IdenticalSequenceWithInSpecies.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("identicalSeq <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(identicalSeq) <- c(\"version\",\"species\", \"molecular\", \"sequence\",\"count\",\"nameIDs\"  )")
	scriptArray.append("## set version factor")
	scriptArray.append("identicalSeq$version <- factor(identicalSeq$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","mature.png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append(" ggplot(identicalSeq[which(identicalSeq$molecular == \"mature\"),],aes(x=version, y=sequence,fill=count)) + geom_tile()+ geom_text(aes(label=count))+ scale_fill_gradient(low=\"yellow\", high=\"red\") +ylab(\"Sequence\")+xlab(\"Version\")+ theme(axis.text.x = element_text(size = 18,hjust = 0,vjust = 0.5,angle=270), legend.position = \"none\" ,axis.title.x = element_text(size=20), axis.text.y = element_text(size = 12),axis.title.y = element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=20), strip.text.x = element_text(size = 18))+facet_wrap(~ species) ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","hairpin.png") +"\",res = 500, units = 'in',width = 30, height = 10)")
	scriptArray.append(" ggplot(identicalSeq[which(identicalSeq$molecular == \"hairpin\"),],aes(x=version, y=sequence,fill=count)) + geom_tile()+ geom_text(aes(label=count))+ scale_fill_gradient(low=\"yellow\", high=\"red\") +ylab(\"Sequence\")+xlab(\"Version\")+ theme(axis.text.x = element_text(size = 18,hjust = 0,vjust = 0.5,angle=270), legend.position = \"none\" ,axis.title.x = element_text(size=20), axis.text.y = element_text(size = 5),axis.title.y = element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=20), strip.text.x = element_text(size = 18))+facet_wrap(~ species) ")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "identicalWithinSpecies.R",  miRscriptsFolder)

def assemblerRscript_NucleotideCombination( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	# inputfile = miOutputFolder+filenamePrefix+"EndNucleotideCombination_mature.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("if (!require('ggseqlogo')) install.packages('ggseqlogo')")
	scriptArray.append("library('ggseqlogo')")
	scriptArray.append("if (!require('scales')) install.packages('scales')")
	scriptArray.append("library('scales')")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("versionFactor <- c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25')")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("## endNucleotide ##")
	scriptArray.append("## read in data")
	inputfile = miOutputFolder+filenamePrefix+"EndNucleotideCombination_mature.tsv"

	scriptArray.append("nucleotideCombination <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(nucleotideCombination) <- c(\"version\",\"species\", \"endPosition\", \"Nucleotides\",\"count\",\"ratio\"  )")
	scriptArray.append("## set version factor")
	scriptArray.append("nucleotideCombination$version <- factor(nucleotideCombination$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("positionlist <- unique(nucleotideCombination$endPosition)")
	scriptArray.append("for (ipos in positionlist){")
	scriptArray.append("iposNTcombination <- nucleotideCombination[which(nucleotideCombination$endPosition == ipos),]")
	scriptArray.append("png(file=gsub(\"tsv\", paste(ipos, \".png\", sep=\"\"),\""+inputfile +"\"),res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append("#")
	scriptArray.append("print( ggplot(iposNTcombination,aes(x=Nucleotides, y=ratio,fill=Nucleotides)) +coord_polar(start = 0)  + geom_bar(stat=\"identity\")+ theme(axis.text.x = element_text(size = 8,hjust = 0,vjust = 0.5,angle=-20), legend.position = \"none\", axis.title.x = element_text(size=20), axis.text.y = element_text(size = 12),axis.title.y = element_text(size=20), strip.text.x = element_text(size = 20)) +facet_grid(version~ species) +ylim(-0.05, max(iposNTcombination$ratio)*1.2) )")
	scriptArray.append("#")
	scriptArray.append("dev.off()")
	scriptArray.append("}")
	scriptArray.append("#")
	scriptArray.append("## seqLogo of seedSeq ##")
	inputfile = miOutputFolder+filenamePrefix+"EndNucleotideCombination_mature.tsv"
	scriptArray.append("## read in data")
	scriptArray.append("seedSeq <- read.table(\""+inputfile+"\",sep = \"\\t\", header = TRUE)")
	scriptArray.append("colnames(seedSeq) <- c(\"version\", \"species\", \"miRNA\", \"seedSeq\" )")
	scriptArray.append("#")
	scriptArray.append("versionlist <- unique(seedSeq$version)")
	scriptArray.append("specieslist <- unique(seedSeq$species)")
	scriptArray.append("#")
	scriptArray.append("for (iversion in versionlist){ ")
	scriptArray.append("filename1 <- paste(\"miRBase_\", iversion,\"_\" ,sep = \"\")")
	scriptArray.append("for (ispecies in specieslist){ ")
	scriptArray.append("filename <- paste(\"" +miOutputFolder+"\", filename1, ispecies, \"_SeedRedionSequence.tsv\", sep = \"\") " )
	scriptArray.append("seedSeq <- read.table(filename,sep = \"\\t\", header = FALSE)")
	scriptArray.append("png(file=gsub(\"tsv\",\"png\",filename),res = 300,units = 'in',width = 10, height = 10)")
	scriptArray.append("print(ggseqlogo(as.character(seedSeq$V2) , method = 'bits', seq_type='rna'))")
	scriptArray.append("dev.off()")
	scriptArray.append("}}")
	scriptArray.append("#")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "NucleotideCombination.R",  miRscriptsFolder)

def assemblerRscript_SeqSimilarityNetwork( ispecies, filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("if (!require('igraph')) install.packages('igraph')")
	scriptArray.append("library('igraph')")

	scriptArray.append("## ")
	scriptArray.append("## read in data: pairwise Levenshtein Distance of miRNAs")
	inputfile = miOutputFolder+filenamePrefix+"mature_"+"LevenshteinDistanceTable.tsv"
	scriptArray.append("matureDistTable <- read.delim(\""+inputfile+"\", sep='\\t',header=FALSE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(matureDistTable) <- c(\"edge1\",\"edge2\", \"dist\"  )")
	scriptArray.append("matureDistTable$color <- lapply(matureDistTable$dist, function(x) colorRampPalette(c(\"darkblue\", \"lightblue\" ))(max(matureDistTable$dist+1))[[x+1]])")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("## read in data: pairwise Levenshtein Distance of pre-miRNAs")
	inputfile = miOutputFolder+filenamePrefix+"hairpin_"+"LevenshteinDistanceTable.tsv"
	scriptArray.append("hairpinDistTable <- read.delim(\""+inputfile+"\", sep='\\t',header=FALSE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(hairpinDistTable) <- c(\"edge1\",\"edge2\", \"dist\"  )")
	scriptArray.append("hairpinDistTable$color <- lapply(hairpinDistTable$dist, function(x) colorRampPalette(c(\"darkgreen\", \"lightgreen\" ))(max(hairpinDistTable$dist+1))[[x+1]])")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("## read in data: chil-parent relationship between miRNAs and pre-miRNAs ")
	inputfile = miOutputFolder+filenamePrefix+"Mature2Hairpin.tsv"
	scriptArray.append("mature2hairpinTable <- read.delim(\""+inputfile+"\", sep='\\t',header=FALSE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(mature2hairpinTable) <- c(\"edge1\",\"edge2\" )")
	scriptArray.append("mature2hairpinTable$dist <- \"-1\"")
	scriptArray.append("mature2hairpinTable$color <- \"red\"")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("edgesNetwork <- rbind(matureDistTable, hairpinDistTable, mature2hairpinTable)")
	scriptArray.append("edgesNetwork$edge1 <- as.character(edgesNetwork$edge1)")
	scriptArray.append("edgesNetwork$edge2 <- as.character(edgesNetwork$edge2)")
	scriptArray.append("edgesNetwork$color <- as.character(edgesNetwork$color)")
	scriptArray.append("#")
	scriptArray.append("## build network graph ")
	scriptArray.append("edgesNetwork.g <- graph.data.frame(edgesNetwork[,c(1,2,4)],directed = FALSE)")
	scriptArray.append("E(edgesNetwork.g)$color <- edgesNetwork$color")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","mature.png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append("#")
	scriptArray.append("plot(edgesNetwork.g, layout= layout.fruchterman.reingold, vertex.size = 1, vertex.color = edgesNetwork$vcolor, vertex.label=NA, vertex.label.cex = 0.5,edge.width=1.5)")
	scriptArray.append("#")
	scriptArray.append("dev.off()")
	scriptArray.append("#")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), ispecies+"-SeqSimilarityNetwork.R",  miRscriptsFolder)

def assemblerRscript_TrackingDiffInVersions( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"IdenticalSequenceWithInSpecies.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("identicalSeq <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(identicalSeq) <- c(\"version\",\"species\", \"molecular\", \"sequence\",\"count\",\"nameIDs\"  )")
	scriptArray.append("## set version factor")
	scriptArray.append("identicalSeq$version <- factor(identicalSeq$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","mature.png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append(" ggplot(identicalSeq[which(identicalSeq$molecular == \"mature\"),],aes(x=version, y=sequence,fill=count)) + geom_tile()+ geom_text(aes(label=count))+ scale_fill_gradient(low=\"yellow\", high=\"red\") +ylab(\"Sequence\")+xlab(\"Version\")+ theme(axis.text.x = element_text(size = 18,hjust = 0,vjust = 0.5,angle=270), legend.position = \"none\" ,axis.title.x = element_text(size=20), axis.text.y = element_text(size = 12),axis.title.y = element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=20), strip.text.x = element_text(size = 18))+facet_wrap(~ species) ")
	scriptArray.append("dev.off()")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "TrackingDiffInVersions.R",  miRscriptsFolder)

def assemblerRscript_HighConfidenceComparison( filenamePrefix, miOutputFolder, miRscriptsFolder):
	scriptArray = []
	inputfile = miOutputFolder+filenamePrefix+"IdenticalSequenceWithInSpecies.tsv"
	scriptArray.append("")
	# scriptArray.append(".libPaths(\"/Users/xiangfuz/Rlibrary/\") ## tmp add")
	scriptArray.append("if (!require('ggplot2')) install.packages('ggplot2')")
	scriptArray.append("library('ggplot2')")
	scriptArray.append("## ")
	scriptArray.append("## read in data")
	scriptArray.append("identicalSeq <- read.delim(\""+inputfile+"\", sep='\\t',header=TRUE)")
	scriptArray.append("## set column name")
	scriptArray.append("colnames(identicalSeq) <- c(\"version\",\"species\", \"molecular\", \"sequence\",\"count\",\"nameIDs\"  )")
	scriptArray.append("## set version factor")
	scriptArray.append("identicalSeq$version <- factor(identicalSeq$version, levels = c( 'V1.0', 'V1.1', 'V1.2', 'V1.3', 'V1.4', 'V1.5', 'V2.0', 'V2.1', 'V2.2', 'V3.0', 'V3.1', 'V4.0', 'V5.0', 'V5.1', 'V6.0', 'V7.0', 'V7.1', 'V8.0', 'V8.1', 'V8.2', 'V9.0', 'V9.1', 'V9.2', 'V10.0', 'V10.1', 'V11.0', 'V12.0', 'V13.0', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', '23', '24', '25'))")
	scriptArray.append("#")
	scriptArray.append("#")
	scriptArray.append("png(file=\""+inputfile.replace("tsv","mature.png") +"\",res = 500, units = 'in',width = 10, height = 10)")
	scriptArray.append(" ggplot(identicalSeq[which(identicalSeq$molecular == \"mature\"),],aes(x=version, y=sequence,fill=count)) + geom_tile()+ geom_text(aes(label=count))+ scale_fill_gradient(low=\"yellow\", high=\"red\") +ylab(\"Sequence\")+xlab(\"Version\")+ theme(axis.text.x = element_text(size = 18,hjust = 0,vjust = 0.5,angle=270), legend.position = \"none\" ,axis.title.x = element_text(size=20), axis.text.y = element_text(size = 12),axis.title.y = element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=20), strip.text.x = element_text(size = 18))+facet_wrap(~ species) ")
	scriptArray.append("dev.off()")
	assemblerRscript_WriteScriptToFile("\n".join(scriptArray), "HighConfidenceComparison.R",  miRscriptsFolder)

def assemblerRscript_WriteScriptToFile( queryList, queryFilename, queryFolder):
	############################################
	## write a list of string, output into file
	############################################
	outputfilename=""
	if queryFolder.endswith("/"):
		outputfilename = queryFolder + queryFilename
	else:
		outputfilename = queryFolder +"/"+ queryFilename
	filehandle = open(outputfilename, "w")
	filehandle.writelines(queryList)
	filehandle.close()
