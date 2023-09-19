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
	#print(SampsD)
	return SampsD



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


def readPersonPhasedFasta(person,personMerge_out_path):
	SeqD,seqLD,RefSeqD = {},{},{}
	##print(personMerge_out_path)
	for File in os.listdir(personMerge_out_path): 
		##print(File)
		if File.split(".")[-1] == "fasta" and "MutatColor" in File:
			Rg = File.split(".fasta")[0]	
			#if "R24453" in Rg :
				##print("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh")		
				##print(Rg)
			SeqD[Rg],seqLD[Rg] = {},{}
			#RgLen = int(Rg.split(".")[2].split("-")[1])
			for Filel in open(SampMerge_out_path + "/" + File).readlines():					
				if ">" in Filel:
					seqID = Filel.split("\n")[0]
				else:
					seq = Filel.split("\n")[0]

					##out the region by region length or the 
					seqLD[Rg] = len(seq)
					#seqLD[Rg] = RgLen
					SeqD[Rg][seqID] = seq
	for Rg in SeqD:
		RefSeqD[Rg] = {}
		for seqid in SeqD[Rg]:
			if ">ref" in seqid:
				RefSeqD[Rg][seqid] = SeqD[Rg][seqid]
	##print(seqLD)
	return SeqD,seqLD,RefSeqD


#Top Rg
def order_dict(dicts, n):
    result = []
    result1 = []
    p = sorted([(k, v) for k, v in dicts.items()], reverse=True)
    s = set()
    for i in p:
        s.add(i[1])
    ##print(p)
    ##print(sorted(s, reverse=True)[:n])
    ##print(len(sorted(s, reverse=True)[:n]))
    ##print("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh")
    for i in sorted(s, reverse=True)[:n]:
        for j in p:
            if j[1] == i:
                result.append(j)
    ct = 0
    for r in result:
        ct += 1
        if ct <= n:
        	result1.append(r[0])
    return result1



def mergeBlock(Rg,SeqD,RefSeqD,maxNumAlloToMerge):	
	##print(Rg)
	seqD = SeqD[Rg]
	#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
		#print(seqD)
		#print("2222222-------------")
	refID = list(RefSeqD[Rg])[0]
	refSeq = RefSeqD[Rg][refID]

	Samp =  person
	SEQsD = {}
	if Samp not in SEQsD:
		SEQsD[Samp] = {}
	sortedSeqs = []
	for seqid in seqD:	
		if "ref_" not in seqid:
			 #seqid.split(">")[1] #.split("-")[0]
			supptNum =  int(seqid.split("\n")[0].split(">")[1].split(".")[0])
			SEQ = seqD[seqid]
			if SEQ not in SEQsD[Samp]:
				SEQsD[Samp][SEQ] = supptNum
			#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
				##print(SEQsD.keys())
				##print(SEQsD["10.3"])
				##print("2222222222222")
			#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
				#print(SEQsD)
				#print("333333333333333----1")

	sorted_sump = sorted(SEQsD[Samp].items(), key=lambda item:item[1], reverse=True)	
	for XXX in range(len(sorted_sump)):
				##print(XXX)
				##print(sorted_sump[XXX])
		sortedSeqs.append(sorted_sump[XXX][0])
	#print(sortedSeqs)
	#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
		#print(sorted_sump)
		#print(sortedSeqs)
		#print("333333333333333")

	##start merge the haps 
	mergeAllSampD,mergeD={},{}
	for seq in sortedSeqs:
		seqid = ''
		for seqDid in seqD:
			if seqD[seqDid] == seq and "ref_" not in seqDid:
				seqid = seqDid
				break
		#print(seqid)
		#print("====")
		if "ref_" not in seqid:
			#Samp =  seqid.split(">")[1] #.split("-")[0]
			supptNum =  int(seqid.split("\n")[0].split(">")[1].split(".")[0])
			SEQ = seqD[seqid]
			
			if Samp not in mergeD:
				mergeD[Samp] = {}

			##print(SEQ)
			if mergeD[Samp] == {}:				
				##print(SEQ)
				Smp_SEQ_sorted_dic = sorted(SEQsD[Samp].items(), key=lambda item:item[1], reverse=True)
				#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
					#print(SEQsD[Samp])
					#print(Smp_SEQ_sorted_dic)
					#print("444444444444---000")
				##print(Smp_SEQ_sorted_dic)
				##print("ooooooooooooo")
				for x in range(1):
					seq = Smp_SEQ_sorted_dic[x][0]
					num = Smp_SEQ_sorted_dic[x][1]
					if num >= 1:
						mergeD[Samp][seq] = num
				#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
					#print(mergeD)
					#print("444444444444")
									
			else:
				#if "R1234923-1247576.L-12654.No-56.MutatColor" in Rg :
					#print(mergeD)
					#print("444444444444----2")
				##print(mergeD)
				##print(SEQ)
				##print("-----------------------------------")
				##print()
				dropSeq,AddSeq,pipeiNo = '','',0
				Count = 0

				Smp_sorted_dic = sorted(mergeD[Samp].items(), key=lambda item:item[1], reverse=True)

				#if "R24453-40974.L-16521.No-384.MutatColor" in Rg :
					##print(Smp_sorted_dic)
					##print("sort")
				##print(Smp_sorted_dic)
				for mgSEQ_K_V in Smp_sorted_dic:
					diffNum,mgSEQBaseNo,SEQBaseNo = 0,0,0
					HebingSeq,HebingSeqBasenum = '',0
					thisSEQ_baseCount = {}
					mgSEQ = mgSEQ_K_V[0]
					##print(mgSEQ)
					for XX in range(len(SEQ)):
						SEQBs = SEQ[XX]
						mgSEQBs = mgSEQ[XX]
						if SEQBs not in thisSEQ_baseCount:
							thisSEQ_baseCount[SEQBs] = 0
						thisSEQ_baseCount[SEQBs] +=1

					SEQ_maxBase = max(thisSEQ_baseCount,key=thisSEQ_baseCount.get)

					for XX in range(len(SEQ)):
						refBs = refSeq[XX]
						SEQBs = SEQ[XX]
						mgSEQBs = mgSEQ[XX]
						if (mgSEQBs != "A" ) and  (SEQBs != "A") and (mgSEQBs != SEQBs):
							diffNum += 1
						if mgSEQBs != "A":
							mgSEQBaseNo += 1							
						if SEQBs != "A" :
							SEQBaseNo += 1

					if len(SEQ) <= 5:
						maxNumAlloToMerge = 0
					else:
						maxNumAlloToMerge = 3
					if len(SEQ) > 0:
						if diffNum <= maxNumAlloToMerge:
							#if mergeD[Samp][mgSEQ] > supptNum:
								#HebingSeq = mgSEQ
							#elif mergeD[Samp][mgSEQ] == supptNum:
							if mgSEQBaseNo >= SEQBaseNo:
								HebingSeq1 = mgSEQ
							else:
								HebingSeq1 = SEQ
							#else:
								#HebingSeq = SEQ
							HebingSeq = ''
							for DX in range(len(HebingSeq1)):
								HebingSeqBase = HebingSeq1[DX]
								Mgbase = mgSEQ[DX]
								ThisSQbase = SEQ[DX]
								if HebingSeqBase == "A" and Mgbase != "A"  :
									HebingSeqBase = Mgbase
								elif HebingSeqBase == "A" and ThisSQbase != "A"  :
									HebingSeqBase = ThisSQbase
								HebingSeq += HebingSeqBase
	#
							#if "R24453-40974.L-16521.No-384.MutatColor" in Rg :
								##print("merge")
								##print(mgSEQ)
								##print(SEQ)
								##print(HebingSeq)
								##print()
								

							dropSeq = mgSEQ
							AddSeq = HebingSeq
							Count += 1
							break
						else:
							#AddSeq = SEQ
							dropSeq = "No"
							#dropSeq = SEQ
					else:
						dropSeq = "No"
						


				if dropSeq != '' and dropSeq != "No":
					##print("drop")
					##print(Count)
					AA = mergeD[Samp][dropSeq]	
					##print(AA)
					##print("hshhs")
					mergeD[Samp].pop(dropSeq)
					mergeD[Samp][AddSeq] = AA
					mergeD[Samp][AddSeq] += supptNum
#
					#if "R24453-40974.L-16521.No-384.MutatColor" in Rg and "BBBBBBBBBBBBBBBBBBBBBBBCCCCCCCCCCCCCCCCBBBBBBBBBBBBBBBBBBBBBC" in SEQ:
						##print("mergege")
						##print(diffNum)
						##print(dropSeq)
						##print(AddSeq)
						##print(mergeD[Samp])
					##print(mergeD)
					##print()
				elif dropSeq == '':
					##print(dropSeq)
					##print("hshhshsh")
					mergeD[Samp][SEQ] = supptNum
				elif dropSeq == 'No':
					mergeD[Samp][SEQ] = supptNum

				else:
					mergeD[Samp][SEQ] += supptNum

							
	######-------------------
	#if(Rg == "R1234923-1247576.L-12654.No-56.MutatColor"):
		#print("mergeD")
		##print(mergeD)
	#if "R24453-" in Rg  and "phasing_P22-C_sL-d" in Samp:
	#	#print(mergeD[Samp])
	##print("ueueusuuuuuuuuuuuuuuuuuuuuuu")
	##print("----------")		
	##print(mergeD)
	##print("----------")
	AllSamp_OK_num = 0			
	for Samp in mergeD:
		FiltHap = {}
		for seq_1 in mergeD[Samp]:
			##print(seq_1)
			##print("=============")
			if mergeD[Samp][seq_1] >= minHapSuptNum:
				FiltHap[seq_1] = mergeD[Samp][seq_1]
				AllSamp_OK_num += 1	
				##print(seq_1)
				##print("hhhhhhhhhhhhhhhhhhhhhhhhhhsssssssssssssssssssss")

		sorted_Seqs = sorted(FiltHap.items(), key=lambda item:item[1], reverse=True)
		#for MergeSEQ in mergeD[Samp]:
			##print("\t".join([Samp,MergeSEQ,str(mergeD[Samp][MergeSEQ]) ] ))
		mergeAllSampD[Samp]	= sorted_Seqs

	##print(mergeAllSampD)
	##print("----------------------0000000000000000000")
	if AllSamp_OK_num == 0 :
		NotOKRg =  "No"
		##print("------------------------")
		##print(Rg)
	else:
		NotOKRg = "yes"


	return NotOKRg,mergeAllSampD


					



def main():
	startTime = time.time()
	################################################################################################################################
	# Check input files
	################################################################################################################################
	SampsFlag = {"t":"t","j":"j","d":"d","x":"x"}
	filtPosiLst,filtPosiD = readPosi(filt_posi_vcfF)
	#print("all Posi Num  :  " + str(len(filtPosiLst)))

	#colorBaseD = {"A":0,"B":1,"C":2,"D":3,"E":4,"F":5,"G":6}
	#colorBaseD = {"A":"A","B":"B","C":"C","D":"D","E":"E","F":"F","G":"G"}
	SampIdDyValue = {"x":5,"t":10,"j":15,"d":20}

	SampsD = readSampInfoF(person)
	##print(SampsD)
	samps  = list(SampsD.keys())
	##print(samps)


	refSamp = person.split("_2_")[1].split("--")[0]
	##print(refSamp)


	RefGnmF = refGnmP + "/" + refSamp + ".fasta"
	RefGnmD  = readRefGnm(RefGnmF)

	RefGnmLst = list(RefGnmD[list(RefGnmD.keys())[0]])

	#print(SampMerge_out_path)
	SeqD,seqLD,RefSeqD = readPersonPhasedFasta(person,SampMerge_out_path)

	#print(SeqD.keys())
	##print(SeqD)
	##print(SeqD["R1234923-1247576.L-12654.No-56.MutatColor"])

	
	AllRegions = seqLD

	AllRegs_merge,OKRegion = {},{}
	for Region in AllRegions:
		##print(Region)
		NotOKRg,mergeAllSampD = mergeBlock(Region,SeqD,RefSeqD,maxNumAlloToMerge)		
		##print(NotOKRg)
		#if(Region == "R1234923-1247576.L-12654.No-56.MutatColor"):
			##print(SeqD[Region])
		
		if NotOKRg != "No":
			##print("NONONONO")
			##print()
			AllRegs_merge[Region] = mergeAllSampD
			OKRegion[Region] = AllRegions[Region]
			##print(mergeAllSampD)
			#if(Region == "R1234923-1247576.L-12654.No-56.MutatColor"):
				##print(mergeAllSampD.keys())	
				##print("---=ahsihsihsi")

	##print(AllRegs_merge)
	Ls = []
	for Rg in AllRegs_merge:
		##print(Rg)
		S = int(Rg.split(".")[0].split("-")[0].split("R")[1])
		E = int(Rg.split(".")[0].split("-")[1])
		L = int(Rg.split(".")[1].split("-")[1])
		
		RegionPosis = []
		for filtPs in filtPosiLst:
			if filtPs >= S and filtPs <= E:
				RegionPosis.append(filtPs)
		##print(RegionPosis)
		
		allNum = 0
		for smp in AllRegs_merge[Rg]:
			##print(smp)			
			for x in range(len(AllRegs_merge[Rg][smp])):	
				#seq = AllRegs_merge[Rg][smp][x][0]
				num = AllRegs_merge[Rg][smp][x][1]
				allNum += num
			##print("========")
			##print(AllRegs_merge[Rg][smp])
			##print(allNum)
			##print("========")
		#if(Rg == "R1234923-1247576.L-12654.No-56.MutatColor"):
			##print("000000=============================")	
			##print(AllRegs_merge[Rg])
			##print(allNum)
		hapNum = 0
		for smp in AllRegs_merge[Rg]:
			#if(Rg == "R1234923-1247576.L-12654.No-56.MutatColor"):
				##print(smp)
				##print(len(AllRegs_merge[Rg][smp]))
				##print("11111111++++")						
			for x in range(len(AllRegs_merge[Rg][smp])):	
				hapNum = hapNum + 1
				#if(Rg == "R1234923-1247576.L-12654.No-56.MutatColor"):
					##print(hapNum)
					##print("qqqqqqqqqqqqqqq")
				seq = AllRegs_merge[Rg][smp][x][0]
				num = AllRegs_merge[Rg][smp][x][1]
				Freq = round(float(num)/allNum,2)

				Start = RegionPosis[0]
				count = 0
				for seqidx in range(len(seq)):
					thisbase = seq[seqidx]
					if seqidx < len(seq)-1:
						nextbase = seq[seqidx +1]
					else:
						nextbase = seq[seqidx]
					if seqidx < len(seq)-2:
						thirdbase = seq[seqidx + 2]
					else:
						thirdbase = seq[seqidx -1]
					
					if seqidx < len(seq) -1:
						if thisbase != nextbase:
							End = RegionPosis[seqidx]								
							base = thisbase							
							outl = []
							outl.append(smp)
							outl.append(Rg)
							outl.append(str(hapNum))
							outl.append(str(num))
							outl.append(str(Freq))
							outl.append(str(Start) + "-" + str(End))
							outl.append(base)
							outl.append(seq)
							#outl.append(str(count))
							l = "\t".join(outl)
							Ls.append(l)

							
							S2 = RegionPosis[seqidx] + 1
							E2 = RegionPosis[seqidx+1] - 1
							#if nextbase == "A":
							#	base2 = "A"	
							#else:
							#	base2 = "B"	
							#if nextbase == "D":
							#	base2 = thisbase	
							#if thisbase == "D":
							#	base2 = nextbase
							#if(thisbase == "A" & (nextbase == "B" | nextbase == "C") ) :
						#		base2 = "A"
						#	if((thisbase == "B" | thisbase == "C") & nextbase == "A"  ) :
						#		base2 = "A"
							if( thirdbase == thisbase) :
								base2 = thisbase
							if( thirdbase == nextbase) :
								base2 = nextbase
							if( thirdbase == nextbase) :
								base2 = nextbase
							if( thirdbase != nextbase and thirdbase != thisbase and thirdbase != "A") :
								base2 = "E"
							if(thisbase == "A"  ) :
								base2 = "A"
							if(nextbase == "A"  ) :
								base2 = "A"
					
								
									
								
							#base2 = "E"	

							outl = []
							outl.append(smp)
							outl.append(Rg)
							outl.append(str(hapNum))
							outl.append(str(num))
							outl.append(str(Freq))
							outl.append(str(S2) + "-" + str(E2))
							outl.append(base2)
							outl.append(seq)
							l = "\t".join(outl)
							Ls.append(l)


							Start = RegionPosis[seqidx+1]
							#count = 0
						#else:
							#count += 1 
					else :
						
						End = RegionPosis[seqidx]
						base = thisbase

						outl = []
						outl.append(smp)
						outl.append(Rg)
						outl.append(str(hapNum))
						outl.append(str(num))
						outl.append(str(Freq))
						outl.append(str(Start) + "-" + str(End))
						outl.append(base)
						outl.append(seq)
						#outl.append(str(count))
						l = "\t".join(outl)
						Ls.append(l)
	##print("\n".join(Ls))
	
	ls = "\n".join(Ls)
	##print(ls)








## out support Nums	of hap line		
	OutHapNumsF = outP + "/" + person + ".mergedHaplotypes.supportNum.txt" 
	if (os.path.exists(OutHapNumsF)) :
		os.remove(OutHapNumsF)
	##print(OutHapNumsF)
	OutHapNumsFO = open(OutHapNumsF,'a')
	OutHapNumsFO.write("\t".join(["#Samp","Region","hapType_id","supportNum","HapFreq","subRegion","Base","haplotype"])+"\n")
	OutHapNumsFO.write("\n".join(Ls)+"\n")
	OutHapNumsFO.close()


## generate fasta file	
	basesDic = {}
	RefGnmLst = list(RefGnmD[list(RefGnmD.keys())[0]])
	#print(RefGnmLst)
	for Cites_AlleBasesFl in open(Cites_AlleBasesF):
		site = int(Cites_AlleBasesFl.split("\t")[1])
		bases = Cites_AlleBasesFl.split("\n")[0].split("\t")[2]
		all_bases = bases.split(";")
		basesDic[site] = {}
		for all_base in all_bases:
			S = all_base.split(":")[0]
			B = all_base.split(":")[1]
			if B == RefGnmLst[int(site)-1]:
				basesDic[site]["ref"] = B
			else:
				basesDic[site][S] = B

	finalhaps,blockHaps = {},{}
	for OutHapNumsl in open(OutHapNumsF).readlines()[1:]:
		Rg = OutHapNumsl.split("\t")[1]
		hap = OutHapNumsl.split("\n")[0].split("\t")[7]
		hap_id =  "hap" +  OutHapNumsl.split("\t")[2] + "_" + OutHapNumsl.split("\t")[3] + "_" + OutHapNumsl.split("\t")[4] + "_" + hap
		if Rg != "Region":	
			S = int(Rg.split(".")[0].split("-")[0].split("R")[1])
			E = int(Rg.split(".")[0].split("-")[1])
			RegionPosis = []
			for filtPs in filtPosiLst:
				if filtPs >= S and filtPs <= E:
					RegionPosis.append(filtPs)
			## ref seq
			refSeq = ''
			for RegionPosi in RegionPosis:
				refSeq += RefGnmLst[int(RegionPosi)-1]

			HAPrealBases = ''
			for SiteIndex in range(0,len(hap)):
				Site = RegionPosis[SiteIndex]
				HapBas = hap[SiteIndex]
				if HapBas == "A":
					realB = "N"
				if HapBas == "B":
					realB = basesDic[Site]["ref"]
				if HapBas == "C":
					realB = basesDic[Site]["2"]
				if HapBas == "D":
					realB = basesDic[Site]["3"]
				if HapBas == "E":
					realB = basesDic[Site]["4"]
				HAPrealBases += realB
					
			finalhap = "\t".join([OutHapNumsl.split("\t")[0],OutHapNumsl.split("\t")[1],\
				OutHapNumsl.split("\t")[2],OutHapNumsl.split("\t")[3],OutHapNumsl.split("\t")[4],\
				hap,HAPrealBases])
			if finalhap not in finalhaps:
				finalhaps[finalhap] = ''
			if Rg not in blockHaps:
				blockHaps[Rg] = {}
				blockHaps[Rg]["ref"] = refSeq
				blockHaps[Rg][hap_id] = HAPrealBases 
			else:
				if OutHapNumsl.split("\t")[2] not in blockHaps[Rg]:
					blockHaps[Rg][hap_id] = HAPrealBases 


	OutHapNumsF = outP + "/" + person + ".mergedHaplotypes.txt" 
	if (os.path.exists(OutHapNumsF)) :
		os.remove(OutHapNumsF)
	##print(OutHapNumsF)
	OutHapNumsFO = open(OutHapNumsF,'a')
	OutHapNumsFO.write("\t".join(["#Samp","Region","hapType_id","supportNum","HapFreq","haplotype_forplot","haplotype"])+"\n")
	OutHapNumsFO.write("\n".join(list(finalhaps.keys()))+"\n")
	OutHapNumsFO.close()


	for RG in blockHaps:
		out = []
		for hap_id in blockHaps[RG]:
			out.append(">" + hap_id)
			out.append(blockHaps[RG][hap_id] )
			
		OutHapNumsF = outP + "/" + ".".join(RG.split(".")[:-1]) + ".mergedHaplotypes.fasta" 
		if (os.path.exists(OutHapNumsF)) :
			os.remove(OutHapNumsF)
		OutHapNumsFO = open(OutHapNumsF,'a')
		OutHapNumsFO.write("\n".join(out)+"\n")
		OutHapNumsFO.close()




if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)


	parser.add_option("-b","--Cites_AlleBasesF",
					  dest = "Cites_AlleBasesF",
					  default = "",
					  metavar = "file",
					  help = "Base composition file of diversity sites [required]")

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
					  
	parser.add_option("-m","--SampMerge_out_path",
					  dest = "SampMerge_out_path",
					  default = "",
					  metavar = "path",
					  help = "merged person's out path  [required]")

	parser.add_option("-o","--outP",
					  dest = "outP",
					  default = "",
					  metavar = "path",
					  help = "out path  [required]")

	parser.add_option("-P","--maxNumAlloToMerge",
					  dest = "maxNumAlloToMerge",
					  default = "10",
					  metavar = "int",
					  help = "max Num Allow the hap To be Merged . [required]")
	parser.add_option("-S","--minHapSuptNum",
					  dest = "minHapSuptNum",
					  default = "5",
					  metavar = "int",
					  help = "min Hap Support Num . [required]")





	(options,args) = parser.parse_args()

	Cites_AlleBasesF	          = os.path.abspath(options.Cites_AlleBasesF)
	filt_posi_vcfF		  = os.path.abspath(options.filt_posi_vcfF)
	SampMerge_out_path        = os.path.abspath(options.SampMerge_out_path)
	refGnmP               = os.path.abspath(options.refGnmP)	
	outP                  = os.path.abspath(options.outP)	
	person                = options.person
	maxNumAlloToMerge = int(options.maxNumAlloToMerge)
	minHapSuptNum = int(options.minHapSuptNum)


	if ( not os.path.exists(outP)):
		os.mkdir(outP)


	main()













