#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time,re 
from optparse import OptionParser
import numpy as np


def readMergeReadSTypes(MergeReadSTypesF):
	SeqD,SeqDep = {},{}
	for MergeReadSTypesFl in open(MergeReadSTypesF).readlines():
		Tags = MergeReadSTypesFl.split("\n")[0].split("\t")
		supNum = int(Tags[0])
		PosiS = Tags[1].split("-")[1]
		PosiE = Tags[1].split("-")[2]
		Seq = ''
		for PoBs in  Tags[2].split("_")[:-1]:
			Po = re.findall("\d+",PoBs)[0]
			#print(Po)
			Bs = PoBs[-1]
			Seq += Bs
			#print(Bs)
		SeqD[PosiS + "_" + PosiE + "_" + Seq] = Seq
		SeqDep[PosiS + "_" + PosiE + "_" + Seq] = supNum


	return SeqD,SeqDep


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



def phasing_search(filtPosiLst,filtPosiD,SeqD,SeqDep,startLociNum,RefGnmLst,refSamp):
	StarSidx = 0
	RightSeq,outPSeq = {},{}
	#for searchIndxS in range(filtPosiD[1546429],len(filtPosiLst)-startLociNum):
	#for searchIndxS in range(15):
	print(len(filtPosiLst)-startLociNum)
	#for searchIndxS in range(len(filtPosiLst)-startLociNum):
	for searchIndxS in range(len(filtPosiLst)-startLociNum):
		if searchIndxS  >= StarSidx:
			print("search start index :" + str(searchIndxS))
			#Ct = 0
			TreeSeqNumD,PhasingFlagD = {},{}
			#for searchIndxE in range(searchIndxS + startLociNum-1,15):
			for searchIndxE in range(searchIndxS+ startLociNum-1,len(filtPosiLst)+1):
			#for searchIndxE in range(searchIndxS + startLociNum,len(filtPosiLst)+1):
				#Ct += 1
				TreeSeqNumD[searchIndxE] = 0
				#print("search end index :" + str(searchIndxE) + "   " + str(filtPosiLst[searchIndxS]) + "-"  + str(filtPosiLst[searchIndxE]))
				searchPosS = filtPosiLst[searchIndxS]
				searchPosE = filtPosiLst[searchIndxE]

				SearchLen  = searchIndxE - searchIndxS +1
				searchPoss = filtPosiLst[searchIndxS:searchIndxE+1]
				if searchPosS == 545833:
					print("SearchLen")
					print(SearchLen)
					print(searchPoss)
					print(searchPosE)
					print()
				TreeSeqNum = 0
				seqBs = {}

				for SeqSE in SeqD:
					RS = int(SeqSE.split("_")[0])
					RE = int(SeqSE.split("_")[1])
					if RS <= searchPoss[0] and RE >= searchPoss[0]:
						RS = searchPoss[0]
					if RE >= searchPoss[-1] and RS <= searchPoss[-1]:
						RE = searchPoss[-1]
						
					RSidx = filtPosiD[RS]
					REidx = filtPosiD[RE]
					SEQ = SeqD[SeqSE]
					Len = len(SeqD[SeqSE])	
					allPosis= filtPosiLst[RSidx:REidx+1]
					if REidx >= searchIndxE:
						CovLen = searchIndxE - RSidx + 1
					else:
						CovLen = REidx - RSidx +1

					#if RS > searchPosE or RE < searchPosS:
					if RS <= searchPosE and RE >= searchPosS:
						if searchPosS == 545833:
							print(RS)
							print(RE)
							print(CovLen)
							print(minCovPerct * SearchLen)
							print(searchPoss)
							print("ooooooooooooooooooooooooo")
							#continue
					if RS <= searchPosE and RE >= searchPosS and CovLen >= minCovPerct * SearchLen:
						if searchIndxS not in RightSeq:
							RightSeq[searchIndxS] = {}
						if searchIndxE not in RightSeq[searchIndxS]:
							RightSeq[searchIndxS][searchIndxE] = {}
						if SeqSE not in RightSeq[searchIndxS][searchIndxE]:
							#print(len(SeqSE))
							RightSeq[searchIndxS][searchIndxE][SeqSE] = SeqDep[SeqSE]

						if searchPosS == 545833:
							print(SEQ)
							print(searchPosE)
							print("shishishishishishihsi")
							print(allPosis)
							print(SearchLen)
							print(len(allPosis))
							#print(RightSeq)
							

						for allPoIndex in range(len(allPosis)):
							allPo = allPosis[allPoIndex]
							Bs = SeqD[SeqSE][allPoIndex]
							if Bs != "-" :
								if allPo not in seqBs:
									seqBs[allPo] = 0
								seqBs[allPo] += 1
						if searchPosS == 545833:
							print(len(seqBs))
							print("len cov")
							print("==========================josjjjjjjjjjjjjj")
							
					#if CovLen >= minCovPerct * SearchLen:
						TreeSeqNumD[searchIndxE] += 1
						#print("seqBs")
						#print(seqBs)

				#print(len(searchPoss))
				#print(seqBs)
				if len(seqBs) == len(searchPoss):
					Flag = "Can"
				else:
					Flag = "Cannot"
				#print(Flag)


				PhasingFlagD[searchIndxE] = Flag
				
				#if searchIndxE -1 >= searchIndxS + startLociNum +1:
				if searchIndxE  >= searchIndxS + startLociNum :

					if searchIndxE-1 in PhasingFlagD and  PhasingFlagD[searchIndxE-1] == "Can" and PhasingFlagD[searchIndxE] == "Cannot": 
					#if PhasingFlagD[searchIndxE] == "Cannot": 
						StarSidx = searchIndxE 
						IDDD = str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS]) + "_" + str(filtPosiLst[searchIndxE-1])
						if IDDD not in outPSeq:
							#print(RightSeq[searchIndxS])
							outPSeq[IDDD] = RightSeq[searchIndxS][searchIndxE-1]
						else:
							outPSeq[IDDD] += RightSeq[searchIndxS][searchIndxE-1]


						if searchPosS == 545833:
							print("111111111111111111111")
							print(Flag)
							print(searchPoss)
							print(seqBs)

						break
					if searchIndxE-1 in PhasingFlagD and  PhasingFlagD[searchIndxE-1] == "Cannot" and PhasingFlagD[searchIndxE] == "Cannot":
						break
					if searchIndxE == len(filtPosiLst)-1 and PhasingFlagD[searchIndxE-1] == "Can" and PhasingFlagD[searchIndxE] == "Can":
						#print("wohahahhahahahha")
						#print(Flag)
						if searchPosS == 545833:
							print("3333333333333333333333333333-==-=-==")
							print(Flag)
						if searchIndxS not in outPSeq:
								
							IDDD = str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS]) + "_" + str(filtPosiLst[searchIndxE])
							if IDDD not in outPSeq:
								outPSeq[IDDD] = RightSeq[searchIndxS][searchIndxE]
							else:
								outPSeq[IDDD] += RightSeq[searchIndxS][searchIndxE]
							

							StarSidx = searchIndxE 
							break

				
					#if startLociNum <= 1 and searchIndxE-1 not in PhasingFlagD and PhasingFlagD[searchIndxE] == "Cannot":
						#StarSidx = searchIndxE 
						#print(searchPosS)
						#if searchPosS == 545833:
							#print("Can - Not")
							#print(searchIndxS)
							##print(RightSeq)
						#IDDD = str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS]) + "_" + str(filtPosiLst[searchIndxE-1])
						#if IDDD not in outPSeq:
							#outPSeq[IDDD] = RightSeq[searchIndxS][searchIndxE-1]
						#else:
							#outPSeq[IDDD] += RightSeq[searchIndxS][searchIndxE-1]
						#break


					if searchIndxE == len(filtPosiLst):
						if searchPosS == 545833:
							print("----------")
							print()
						IDDD = str(searchIndxS) + "_" + str(filtPosiLst[searchIndxS]) + "_" + str(filtPosiLst[searchIndxE])
						if IDDD not in outPSeq:
							outPSeq[IDDD] = RightSeq[searchIndxS][searchIndxE]
						else:
							outPSeq[IDDD] += RightSeq[searchIndxS][searchIndxE]
						break

										
	#print(outPSeq)
	#print()
	#print("ooooooooooooooooooooooooooo")
	Lenslst,containPosiNum,Outsummaryls = [],0,[]
	for RegionID in outPSeq:
		#endPo = list(outPSeq[RegionID].keys())[0].split("_")[1]
		startPo = RegionID.split("_")[1]
		endPo = RegionID.split("_")[2]
		Sidx = filtPosiD[int(startPo)]
		Eidx = filtPosiD[int(endPo)]
		#print(startPo)
		#print(endPo)
		
		RgAllPosis = filtPosiLst[Sidx:Eidx+1]
		refSEQ = ''
		for refP in RgAllPosis:
			refBase = RefGnmLst[int(refP)-1]
			refSEQ += refBase
		

		OUTD,validOUTCout,outID_Dic,depthL,TreeSeqNum,SeqLstForSort,outseqD = {},0,{},[],0,[],{}
		for RegionSeqID in outPSeq[RegionID]:
			rgRdS = RegionSeqID.split("_")[0]
			rgRdE = RegionSeqID.split("_")[1]
			RdSeq = RegionSeqID.split("_")[2]
			rgRdSidx = filtPosiD[int(rgRdS)]
			rgRdEidx = filtPosiD[int(rgRdE)]
			rgRdPosiNum = rgRdEidx - rgRdSidx +1
			Depth = outPSeq[RegionID][RegionSeqID]
			RgRdAllPosis = filtPosiLst[rgRdSidx:rgRdEidx+1]
			
			RgRdAllPoBases = {}
			for RgRdAllPoIdx in range(len(RgRdAllPosis)):
				base = RdSeq[RgRdAllPoIdx]
				RgRdAllPo = RgRdAllPosis[RgRdAllPoIdx]
				RgRdAllPoBases[RgRdAllPo] = base

			#print(RgRdAllPoBases)
			outseq,BseCt = '',0
			for RgAllPosi in RgAllPosis:
				if RgAllPosi in RgRdAllPosis:
					Base = RgRdAllPoBases[RgAllPosi]
					BseCt += 1
				else:
					Base = 'N'
				outseq += Base
			#outID = "R" + rgRdS + "-" + rgRdE + ".L_" + str(int(rgRdE)-int(rgRdS))  + ".PosNo_" + str(rgRdPosiNum) + ".D_" + Depth
			#outID = SampsFlag[Sample.split("-")[-1]] + "-" + Depth 
			if BseCt >= len(RgAllPosis) * minCovPerct:
				Flg = outseq #+ "__" + SampsFlag[Sample.split("-")[-1]]
				if Flg not in outseqD:
					outseqD[ Flg ] = Depth
				else:
					outseqD[ Flg ] += Depth



		for SEQ in outseqD:
			SEQdepth = outseqD[SEQ]
			if int(SEQdepth) >= minDepthOut :
				outID = str(SEQdepth) + "." + str(validOUTCout)
				#outID = SampsFlag[Sample.split("-")[-1]] + "-" + str(SEQdepth) + "." + str(validOUTCout)

				SEQout = SEQ.split("__")[0]
				validOUTCout += 1
				OUTD[outID ] =  SEQout  #+ "-" + str(outID_Dic[outID])				
				depthL.append(int(SEQdepth))
				#print(SEQout)
				SeqLstForSort.append(SEQout + ".." + outID ) #+ "-" + str(outID_Dic[outID])

		#print(SeqLstForSort)
		#print(OUTD)
		if validOUTCout >= minhapNum_forOut:
			Ols = []

			Ols.append(">ref_" + refSamp )    
			Ols.append(refSEQ)
			SeqLstForSort.sort()
			for IDidx in range(len(SeqLstForSort)):
				ID = SeqLstForSort[IDidx].split("..")[1]
				Ols.append(">" + ID )    #+ "." + str(IDidx) 
				Ols.append(OUTD[ID])
				
				
			PosiNum = Eidx - Sidx +1
			containPosiNum += PosiNum
			Lenslst.append(int(endPo)-int(startPo))
			Outsummaryls.append("\t".join([startPo,endPo,str(int(endPo)-int(startPo)+1),str(PosiNum),str(np.sum(depthL))]))

			outF = outP + "/R" + startPo + "-" + endPo  + ".L-" + str(int(endPo)-int(startPo)+1) + ".No-" + str(PosiNum) + ".fasta"

			O = "\n".join(Ols)


			## Full info in fasta
			if (os.path.exists(outF)) :
				os.remove(outF)
			ONT_FO = open(outF,'a')
			ONT_FO.write(O+"\n")
			ONT_FO.close()




##summary
		StatF = outP + "/Phasing_summary.txt"
		Len_mean = round(np.mean(Lenslst),2)
		Len_sum = np.sum(Lenslst)
		Lenslst.sort()
		Len_median = np.median(Lenslst)
		Len_max = np.max(Lenslst)

		Samp = outP.split("/")[-1].split("phasing_")[1].split("_2_")[0]
		RefSamp = outP.split("/")[-1].split("phasing_")[1].split("_2_")[1].split("--")[0]
		minFreq = outP.split("/")[-1].split("--")[1].split("-")[0]
		maxFreq = outP.split("/")[-1].split("--")[1].split("-")[1]


		summary="\t".join([Samp,RefSamp,minFreq,maxFreq,str(len(Lenslst)),str(Len_mean),str(Len_median),str(Len_max),str(Len_sum),str(containPosiNum),str(len(filtPosiLst)),str(float(round(containPosiNum *100/len(filtPosiLst),2)))])
		Head = "\t".join(["Samp","RefSamp","FreqMin","FreqMax","contigNum","MeanLen","MedianLen","maxLen","Len_sum","PhasedLociNum","allReliaLociNum","percentLociNum(%)"])
		if (os.path.exists(StatF)) :
			os.remove(StatF)
		ONT_FO = open(StatF,'a')
		ONT_FO.write(Head+"\n")
		ONT_FO.write(summary+"\n")
		ONT_FO.write("\n")
		ONT_FO.write("\n")

		ONT_FO.write("Phased Region summary : "+ "\n")
		ONT_FO.write("\n")
		
		RegionSummaryH="\t".join(["StartLoci","EndLoci","RegionLen","PhasedLociNum","DepthSum"])
		RegionSummaryls = "\n".join(Outsummaryls)

		ONT_FO.write(RegionSummaryH+"\n")
		ONT_FO.write(RegionSummaryls+"\n")

		ONT_FO.close()






def main():
	startTime = time.time()
	################################################################################################################################
	# Check input files
	################################################################################################################################
	#SampsFlag = {"t":"t","j":"j","d":"d","x":"x"}
	filtPosiLst,filtPosiD = readPosi(filt_posi_vcfF)
	print("all Posi Num  :  " + str(len(filtPosiLst)))
	print(MergeReadSTypesF)
	SeqD,SeqDep = readMergeReadSTypes(MergeReadSTypesF)
	#print(SeqD)

	RefGnmD  = readRefGnm(refSampGnmF)
	RefGnmLst = list(RefGnmD[list(RefGnmD.keys())[0]])


	phasing_search(filtPosiLst,filtPosiD,SeqD,SeqDep,startLociNum,RefGnmLst,refSamp)



##log 
	if (os.path.exists(out_line_file)) :
		os.remove(out_line_file)
	ONT_FO = open(out_line_file,'a')
	ONT_FO.write("step:3"+"\n")
	ONT_FO.close()


	endTime = time.time()
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")








if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)


	parser.add_option("-R","--refSampGnmF",
					  dest = "refSampGnmF",
					  default = "",
					  metavar = "file",
					  help = "refSamp Gnm file [required]")
					  
	parser.add_option("-f","--filt_posi_vcfF",
					  dest = "filt_posi_vcfF",
					  default = "",
					  metavar = "string",
					  help = "filt posi File [required]")

	parser.add_option("-F","--MergeReadSTypesF",
					  dest = "MergeReadSTypesF",
					  default = "",
					  metavar = "file",
					  help = "out path  [required]")

	parser.add_option("-o","--outP",
					  dest = "outP",
					  default = "",
					  metavar = "path",
					  help = "out path  [required]")

	parser.add_option("-l","--out_line_file",
					  dest = "out_line_file",
					  default = "",
					  metavar = "file",
					  help = "out_line_file . [required]")
	parser.add_option("-s","--startLociNum",
					  dest = "startLociNum",
					  default = "",
					  metavar = "int",
					  help = "start search Loci Num . [required]")
	parser.add_option("-m","--minhapNum_forOut",
					  dest = "minhapNum_forOut",
					  default = "",
					  metavar = "int",
					  help = "min Loci Num for Output . [required]")
	parser.add_option("-D","--minDepthOut",
					  dest = "minDepthOut",
					  default = "",
					  metavar = "int",
					  help = "min Depth of haplotype to Output . [required]")
	parser.add_option("-c","--minCovPerct",
					  dest = "minCovPerct",
					  default = "0.9",
					  metavar = "float",
					  help = "min num of base num of the read in the phasing region  . [required]")

	parser.add_option("-S","--Sample",
					  dest = "Sample",
					  default = "",
					  metavar = "str",
					  help = "Sample . [required]")
	parser.add_option("-r","--refSamp",
					  dest = "refSamp",
					  default = "",
					  metavar = "str",
					  help = "refSamp . [required]")

	(options,args) = parser.parse_args()

	refSampGnmF	          = os.path.abspath(options.refSampGnmF)
	filt_posi_vcfF		  = os.path.abspath(options.filt_posi_vcfF)
	MergeReadSTypesF      = os.path.abspath(options.MergeReadSTypesF)	
	outP                  = os.path.abspath(options.outP)	
	out_line_file	      = os.path.abspath(options.out_line_file)
	startLociNum	      = int(options.startLociNum)
	minhapNum_forOut	  = int(options.minhapNum_forOut)
	minCovPerct           = float(options.minCovPerct)
	Sample                = options.Sample
	refSamp                = options.refSamp
	minDepthOut           = int(options.minDepthOut)

	main()













