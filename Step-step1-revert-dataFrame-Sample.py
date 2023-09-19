#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time,re 
from optparse import OptionParser
import numpy as np


def readSampInfoF(person):
	SampsD = {}
	Samp = person.split("_2_")[0]
	refSamp = person.split("_2_")[1]
	SampsD[Samp] = refSamp
	print(SampsD)
	return SampsD


def readSampInfoF1(person,SampInfoF):
	SampsD,SampFreq = {},{}
	for SampInfoFl in open(SampInfoF).readlines():
		Tags = SampInfoFl.split("\n")[0].split("\t")
		Samp = Tags[0]
		refSamp = Tags[1]
		personID = Tags[0]
		SampFreqs = Tags[3] + "-" + Tags[4]

		print(personID)
		print(person)
		print(person.split("phasing_")[1].split("_2_")[0])
		print()

		if personID == person.split("phasing_")[1].split("_2_")[0]:
			print("hshhshs")
			SampsD[Samp] = refSamp
			SampFreq[Samp] =SampFreqs
	print(SampsD)
	return SampsD,SampFreq


def readPosi(filt_posi_vcfF):	
	filtPosiLst,filtPosiD = [],{}
	filt_posi_ls = open(filt_posi_vcfF,'r').readlines()
	Count = 0
	for filt_posi_l in filt_posi_ls:
		if "#" not in filt_posi_l:
			filtPosiLst.append(int(filt_posi_l.split("\t")[0]))
			filtPosiD[int(filt_posi_l.split("\t")[0])] = Count
			Count += 1
	return filtPosiLst,filtPosiD



def readRefGnm(RefGnmF):
	RefGnm = {}
	for RefGnmFl in open(RefGnmF).readlines():
		if ">" in RefGnmFl:
			ID = RefGnmFl.split("\n")[0]
			seq = ''
		else:
			seq += RefGnmFl.split("\n")[0]
	RefGnm[ID] = seq
	return RefGnm




def readPersonPhasedFasta(SampMerge_out_path):
	SeqD,seqLD,Nindex,uniqSeqD = {},{},{},{}
	for File in os.listdir(SampMerge_out_path): 
		if File.split(".")[-1] == "fasta" and ".MutatColor" not in File:
			Rg = File.split(".fasta")[0]
			SeqD[Rg],seqLD[Rg] = {},{}

			for Filel in open(SampMerge_out_path + "/" + File).readlines():					
				if ">" in Filel:
					seqID = Filel.split("\n")[0]
				else:
					seq = Filel.split("\n")[0]
					seqLD[Rg] = len(seq)
					SeqD[Rg][seqID] = seq
	
	for Rg in SeqD:
		Nindex[Rg] = {}
		for idx in range(seqLD[Rg]):
			Nindex[Rg][idx] = 0
			for seqid in SeqD[Rg]:
				if SeqD[Rg][seqid][idx] == "N" : # or SeqD[Rg][seqid][idx] == "-":
					Nindex[Rg][idx] += 1


	return SeqD,Nindex





def main():
	startTime = time.time()
	################################################################################################################################
	# Check input files
	################################################################################################################################
	
	filtPosiLst,filtPosiD = readPosi(filt_posi_vcfF)
	print("all Posi Num  :  " + str(len(filtPosiLst)))


	SampsD = readSampInfoF(person)
	samps  = list(SampsD.keys())

	print(samps)
	refSamp = SampsD[samps[0]]
	RefGnmF = refGnmP + "/" + refSamp + ".fasta"
	RefGnmD  = readRefGnm(RefGnmF)
	RefGnmLst = list(RefGnmD[list(RefGnmD.keys())[0]])
	#print(RefGnmLst)
	print(merge_out_path)
	SampMerge_out_path = merge_out_path
	SeqD,Nindex = readPersonPhasedFasta(SampMerge_out_path)

	print(SeqD.keys())
	print(SampMerge_out_path)
	colorBaseD = {0:"A",1:"B",2:"C",3:"D",4:"E",5:"F",6:"G",8:"X"}
	
	SampColorD = {}
	for Region in SeqD:
		RS = int(Region.split(".")[0].split("-")[0].split("R")[1])
		RE = int(Region.split(".")[0].split("-")[1])
		RegionPoss = []
		for FiltPosi in filtPosiLst:
			if FiltPosi >= RS and FiltPosi <= RE:
				RegionPoss.append(FiltPosi)
		RefSeq = ''
		for RegionPos in RegionPoss:
			RefBs = RefGnmLst[RegionPos-1]
			RefSeq += RefBs
		print("ref-REF")
		print(RefSeq)
		num = 0
		compSeqPosiBases = {}
		for RefSeqidx in range(len(RefSeq)):
			#Posi = RegionPoss[RefSeqidx]
			RefBase = RefSeq[RefSeqidx]
			compSeqPosiBases[RefSeqidx] = {1:RefBase}

		OutRegionMutColor,OutSampMutColor,OutPosiDist,OutGroupInfo = [],[],[],[]
		#OutRegionMutColor.append(">Ref_" + refSamp)		
		#OutRegionMutColor.append(RefSeq)	
		OutSampMutColor.append( "\t".join(["SEQid","PosIdx","MutatColor"]) )
		OutPosiDist.append("\t".join(["posIdx","distance"]))
		OutGroupInfo.append("\t".join(["SEQid","group"]))
		for PosIdx in range(len(RegionPoss)-1):
			nowPos = RegionPoss[PosIdx]
			nextPos = RegionPoss[PosIdx+1]
			dist = int(nextPos) - int(nowPos)
			OutPosiDist.append("\t".join([str(PosIdx +1),str(dist)]))




		for seqID in SeqD[Region]:
			#print(seqID)
			if num <= len(SeqD[Region]):
				seq = SeqD[Region][seqID]
				MutateColorLst,Outl = [],[]
				#print(len(seq))
				#OutSampMutColor.append( "\t".join([seqID.split(">")[1],"1","8"]) )
				for seqidx in range(len(seq)):
					seqBase = seq[seqidx]
					#RefBase = RefSeq[seqidx]
					sampColorl = [seqID.split(">")[1]]
					#if seqBase == ".":
						#color_mutat = 0
					#if "ref_" in seqID:
						#color_mutat = 8
					
					if seqBase == "N" or seqBase == "-":
						color_mutat = 0
					else:
						findFlag = "no"
						for CompBsid in compSeqPosiBases[seqidx]:
							if seqBase == compSeqPosiBases[seqidx][CompBsid] :
								color_mutat = CompBsid 
								findFlag = "yes"
								break
						if findFlag == "no":
							color_mutat = len(compSeqPosiBases[seqidx]) + 1
							compSeqPosiBases[seqidx][len(compSeqPosiBases[seqidx]) + 1] = seqBase
					
					#print(compSeqPosiBases)
					shijiPosi = RegionPoss[seqidx]
					if shijiPosi not in SampColorD:
						SampColorD[shijiPosi] = {}
					if seqBase != "N" and seqBase != "-":
						SampColorD[shijiPosi][color_mutat] = seqBase
					#print(shijiPosi)
					#print(color_mutat)
					#print()
					MutateColorLst.append(colorBaseD[color_mutat])
					#MutateColorLst.append(str(color_mutat))
					sampColorl.append(str(seqidx+1))
					sampColorl.append(str(color_mutat))
					OutSampMutColor.append("\t".join(sampColorl))	
					
					if ">x-" in seqID and shijiPosi == 2801:
						print(SampColorD[shijiPosi])
						print("oooooooooooooooooooooo")
				#for poidx in range(len(MutateColorLst)):
					##print(colorBaseD)
					#ColNum = MutateColorLst[poidx]
					##print(ColNum)
					#ColoBase = colorBaseD[ColNum]
					#Outl.append(ColoBase)
				
				OutRegionMutColor.append(seqID)		
				OutRegionMutColor.append("".join(MutateColorLst))

			OutGroupInfo.append("\t".join([seqID.split(">")[1],seqID.split(">")[1].split("-")[0]]))		
				#OutRegionMutColor.append("0" + "".join(MutateColorLst) + "0")		
			num += 1


		Out = "\n".join(OutRegionMutColor)	
			
		outF = outP + "/" + Region + ".MutatColor.fasta" 
		if (os.path.exists(outF)) :
			os.remove(outF)
		#print(outF)
		ONT_FO = open(outF,'a')
		ONT_FO.write(Out+"\n")
		ONT_FO.close()		


		# Out = "\n".join(OutGroupInfo)	
		# outF = outP + "/" + Region + ".groupInfo.txt" 
		# if (os.path.exists(outF)) :
		# 	os.remove(outF)
		# #print(outF)
		# ONT_FO = open(outF,'a')
		# ONT_FO.write(Out+"\n")
		# ONT_FO.close()		


		# OutC = "\n".join(OutSampMutColor)				
		# outF = outP + "/" + Region + ".MutatColor.txt" 
		# if (os.path.exists(outF)) :
		# 	os.remove(outF)
		# #print(outF)
		# ONT_FO = open(outF,'a')
		# ONT_FO.write(OutC+"\n")
		# ONT_FO.close()		


		# OutDist = "\n".join(OutPosiDist)				
		# outDF = outP + "/" + Region + ".distance.txt" 
		# if (os.path.exists(outDF)) :
		# 	os.remove(outDF)
		# #print(outDF)
		# ONT_FO = open(outDF,'a')
		# ONT_FO.write(OutDist+"\n")
		# ONT_FO.close()		



	## out  sample bases in each position
	outF = outP + "/" + samps[0].split("_2_")[0] + ".Posi_AlleBases.txt" 
	if (os.path.exists(outF)) :
		os.remove(outF)
	#print(outF)

	#print(SampColorD)
	lsL = []
	#if ">x-" in seqID and shijiPosi == 2801:
	#print(SampColorD[2801])
	for shijiPosi in SampColorD:
		lL = [samps[0].split("_2_")[0],str(shijiPosi)]
		lTag = []
		for alleOrder in range(1,len(SampColorD[shijiPosi])+1 ):
			if alleOrder in SampColorD[shijiPosi]:
				lTag.append(str(alleOrder) + ":" + SampColorD[shijiPosi][alleOrder])
		lL.append(";".join(lTag))
		if lTag != []:
			lsL.append("\t".join(lL))
	Os = "\n".join(lsL)
	ONT_FO = open(outF,'a')
	ONT_FO.write(Os+"\n")
	ONT_FO.close()		






if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)


	# parser.add_option("-I","--SampInfoF",
	# 				  dest = "SampInfoF",
	# 				  default = "",
	# 				  metavar = "file",
	# 				  help = "Samp  Config info file [required]")

	parser.add_option("-v","--filt_posi_vcfF",
					  dest = "filt_posi_vcfF",
					  default = "",
					  metavar = "string",
					  help = "filt posi File [required]")

	parser.add_option("-r","--refGnmP",
					  dest = "refGnmP",
					  default = "",
					  metavar = "path",
					  help = "ref Gnm P [required]")

	parser.add_option("-p","--person",
					  dest = "person",
					  default = "",
					  metavar = "str",
					  help = "person ID [required]")
					  
	parser.add_option("-m","--merge_out_path",
					  dest = "merge_out_path",
					  default = "",
					  metavar = "path",
					  help = "merged person's out path  [required]")

	parser.add_option("-o","--outP",
					  dest = "outP",
					  default = "",
					  metavar = "path",
					  help = "out path  [required]")


	#parser.add_option("-m","--minDepth",
					  #dest = "minDepth",
					  #default = "1",
					  #metavar = "int",
					  #help = "min Depth to output the haplotype . [required]")


	(options,args) = parser.parse_args()

	#config_F	          = os.path.abspath(options.config_F)
	filt_posi_vcfF		  = os.path.abspath(options.filt_posi_vcfF)
	#SampInfoF             = os.path.abspath(options.SampInfoF)	
	merge_out_path        = os.path.abspath(options.merge_out_path)
	refGnmP               = os.path.abspath(options.refGnmP)	
	outP                  = os.path.abspath(options.outP)	
	person                = options.person

	if ( not os.path.exists(outP)):
		os.mkdir(outP)


	main()













