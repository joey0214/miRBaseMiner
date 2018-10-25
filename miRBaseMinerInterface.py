#! /usr/bin/python

## programmed in Python 2.7

## Author: Xiangfu (Joey) Zhong
## Any feedback can send to joey.zhong.cn@gmail.com

import os
import miRBaseMiner as miner

####################################################
####################################################
## DEMOSTRATION ## 
####################################################
# miner.Version()
# miner.DMOminning()

####################################################
####################################################
## SET UP ##
####################################################
## Please set work directory for miRBaseMiner. If not specifiy, it will work on current work directory
miRBaseMining_workDir = "/PATH_TO"

miRBaseMining_targetVersion = "22" ## ONLY single version

miRBaseMining_versionList = ["9.2","10.1","11.0","12.0","13.0","14","15","16","17","18","19","20","21","22"]

miRBaseMining_species = ["hsa","mmu"] ## can also be list ["hsa", "mmu"]

miRBaseMining_outputFolder = "" ## if not specifiy, it will create a new folder "resultsTables" in current work dir

miRBaseMining_NTmers = ["2","3"]  ## set it as list for furtuer extension 

miRBaseMining_LevenshteinDistanceCutoff = 3 ## the default is set to 999, if want to specifify it, set it to the value you want, (e.g. 3), result will list miRNAs that their distance is less than the cutoff. check any miRNA <= cutoff

## option for building a curated set of miRNA annotation
# miRBaseMining_gffOptions = ["RCm", "RCp", "oRCm", "oRCp", "3LDmp", "ISm", "OSm", "ISp", "OSp","SEQUENCE","NAME", "HF"]
miRBaseMining_gffOptions = ["3LDmp", "ISm", "ISp", "OSm", "OSp"]
## key word for build curated annotation
## RCm	---	full length reverse complementary in mature
## RCp	---	full length reverse complementary in precursor
## oRCm	---	Overlapping reverse complementary in mature
## oRCp	---	Overlapping reverse complementary in precursor
## 0LD	---	maximum allowed pairwise Levenshtein Distance: 0. m p mp
## 1LD	---	maximum allowed pairwise Levenshtein Distance: 1. m p mp
## 2LD	---	maximum allowed pairwise Levenshtein Distance: 2. m p mp
## 3LD	---	maximum allowed pairwise Levenshtein Distance: 3. m p mp
## ISm	---	Identical sequence in precursor sequence
## OSm	---	Overlapping sequence in precursor 
## ISp	---	Identical sequence in mature sequence
## OSp	---	Overlapping sequence in precursor sequence
## -20MFE-80	---	precursor MFE range
## 19ML22	---	mature sequence length range
## 50PL70	---	precursor sequence length range
## GCOR	---	genome coordinate
## HF	---	only including high confidence annotation
## NEW	---	entries that newly added
## NAME	---	entries that changed name
## SEQUENCE	---	entries that changed sequence

####################################################
####################################################
## parameters checking ##
####################################################
if miRBaseMining_workDir == "":
	miRBaseMining_workDir = os.getcwd()
elif not os.path.exists(miRBaseMining_workDir):
	os.mkdir(miRBaseMining_workDir)
miRBaseMining_LevenshteinDistanceCutoff = int(miRBaseMining_LevenshteinDistanceCutoff) # in case, user put a string here



####################################################
####################################################
## main program ##
####################################################


# miner.MFEdensityInMiRStr(miRBaseMining_workDir, ["20", "21", "22"], miRBaseMining_outputFolder)

#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate
#### OUTPUT: all accessible files within each chosen miRBase release
#### DESC: Download data from miRBase ftp site. Dependency: (1) wget; (2) pigz. 
# miner.RetrieveMiRBase(miRBaseMining_workDir, miRBaseMining_versionList)

#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder
#### OUTPUT: (1) the amount of miRNAs and pre-miRNAs in each species in every version
#### OUTPUT: (2) sequence length distribution of mature miRNAs
#### OUTPUT: (3) the varation between amount of miRNA and pre-miRNA in each species
#### OUTPUT: (4) and word cloud for species abundancy in each version
#### DESC:  investigate the overview of miRBase 

#miner.Overview(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder)


#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT: key word in output file name "fullDistribution4MFE" 
#### DESC: investigate MFE distribution of pre-miRNA. also including other features related to MFE: sequence length, GC content, MFE index
# miner.MFEdistribution(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species)

#### identical sequence miRNAs
#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT: (1) profile about identical sequence in miRNAs in each sspecies, key "denticalSequenceBetweenSpecies"; (2) profile about miRNA sequence share between species, "IdenticalSequenceWithInSpecies"
#### DESC: (1) investigate full length identical sequence in miRNA/pre-miRNA, share sequence between different species. require more than one species in miSpecies. 
#### DESC: (2) BOTH for miRNA and pre-miRNA
# miner.IdenticalSequenceMir(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species)


#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT: "SequenceOverlapping" 
#### DESC:  find sequence overlapping with others in given species list and version list. In miRNA and pre-miRNA
# miner.SequenceOverlapping(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species)

#### investigate Nucleotide combination
#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species; (5) the nucleotide mer number.
#### OUTPUT:
#### DESC:  investigate the nucleotide combination of miRNA. ONLY available for mature miRNAs. 
#### DESC: (1) endNucleotide; (2) polyAtailing; (3) polyNTheading; (4) seedRegionSeq.
# miner.NucleotideCombination(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species, miRBaseMining_NTmers)


#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT:
#### DESC:  tracking what happen in "miRNA.diff" file; can catch some entries that have multiple change tags (i.e. NAME and SEQUENCE)  
#### DESC: summarize :
#### DESC:  (1)how many NEW miRNAs were added  
#### DESC:  (2)how many miRNAs were  DELETE 
#### DESC:  (3)how many miRNAs were  changed SEQUENCE
#### DESC:  (4)how many miRNAs were  changed NAME
#### DESC:  NEW ---------- new to the database
#### DESC:  DELETE ------- removed since the last release
#### DESC:  SEQUENCE ----- sequence has changed
#### DESC:  NAME --------- ID has changed
# miner.TrackingDiffInVersions(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species )

 
#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT:
#### DESC: investigate ReverseComplementary. reverse complementary within each given species in given versions. only focus on the sequence, regardless the sequence annotated on which strand
# miner.ReverseComplementary(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species )

 
#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species.
#### OUTPUT:
#### DESC:  muliple locations for miRNAs and pre-miRNAs. this only works if the version has subfolder (genome) and a gff file
#### DESC:  collect genome coordinates for miRNA/pre-miRNA from gff annotation file
#### DESC:  use miRNA/pre-miRNA ID/Accession as key to pull multiple location information together
# miner.MiRrMultipleLocations(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species )


#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species; (5) the target version of miRBase release want to compre.
#### OUTPUT:
#### DESC: compare to high confidence annotation set
#### DESC: this functon is hard coded, currently work on release v20, v21, v22. if you wish to work on later version (if available ), then you nedd to access to MinerUtility.py to specify the file name it should work on that version
# miner.HighConfidenceComparison(miRBaseMining_workDir, miRBaseMining_versionList, miRBaseMining_outputFolder, miRBaseMining_species, miRBaseMining_targetVersion )


#### INPUT: (1) the work directory for miRBaseMiner; (2) versions that want to investigate; (3) output folder; (4) interesting species; (5) the cutoff for pairwise Levenshtein distance.
#### OUTPUT:
#### DESC: similarity network. miRBaseMining_LevenshteinDistanceCutoff will apply to both miRNA and pre-miRNA
# miner.SeqSimilarityNetwork(miRBaseMining_workDir, miRBaseMining_targetVersion, miRBaseMining_outputFolder, miRBaseMining_species, miRBaseMining_LevenshteinDistanceCutoff)


#### INPUT: (1) the work directory for miRBaseMiner; (2) the base release, which curated annotation based on; (3) output folder; (4) interesting species; (5) versions that want to investigate;  (6) options for building curated annotation set.
#### OUTPUT:
#### DESC: Build curated miRNA annotation
#### DESC: the miRBaseMining_versionList is indendpend from other function. The versions that you want to coonsider for building curated annotation profile
#### DESC: the entries without genome coordinate will automatically exclude in curated gff file 
# miner.BuildingCuratedGff(miRBaseMining_workDir, miRBaseMining_targetVersion, miRBaseMining_outputFolder, miRBaseMining_species, miRBaseMining_versionList, miRBaseMining_gffOptions)
############ keywords curated annotation options ############
## key word for build curated annotation
## RCm	---	full length reverse complementary in mature
## RCp	---	full length reverse complementary in precursor
## oRCm	---	Overlapping reverse complementary in mature
## oRCp	---	Overlapping reverse complementary in precursor
## 0LD	---	maximum allowed pairwise Levenshtein Distance: 0. m p mp
## 1LD	---	maximum allowed pairwise Levenshtein Distance: 1. m p mp
## 2LD	---	maximum allowed pairwise Levenshtein Distance: 2. m p mp
## 3LD	---	maximum allowed pairwise Levenshtein Distance: 3. m p mp
## ISm	---	Identical sequence in precursor sequence
## OSm	---	Overlapping sequence in precursor 
## ISp	---	Identical sequence in mature sequence
## OSp	---	Overlapping sequence in precursor sequence
## -20MFE-80	---	precursor MFE range
## 19ML22	---	mature sequence length range
## 50PL70	---	precursor sequence length range
## GCOR	---	genome coordinate
## HF	---	only including high confidence annotation
## NEW
## NAME
## SEQUENCE



