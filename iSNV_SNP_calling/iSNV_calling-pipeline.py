#coding=UTF-8
###################################
#######iSNV_calling 20210517
#######
###################################



personLst = ["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25"]

#personLst = ["P1","P2"]



rule all:
	input:
	#bwa map
		## iSNV calling
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/all.iSNV_with_SNP.pyResults.txt",person=personLst),
		##vcf
		expand("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/vcf",person=personLst),


rule ntfreq_2_FreqBigtable:
	input:
		configF="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/table_config/{person}.ntfreq_file_list.txt",  ##format in the config file: path/sampID.ntfreq"\n"
		refF="/shared/liuhj/HP/process/assembly/person_respectGnm/{person}.respect.fasta",
	output:	
		outputF="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		outputF2="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/out.samps-ntfreq.bigtable.txt",

	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		outputP="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}",
		DEP_THRES=100,       ##最低深度
		MIN_DEP_THRES=5,    ##支持alt的最小reads数
		VALID_SIZE_THRES=1000000,  ##样本有效需要的有效位点数
		FREQ_THRES=0.02,    ###alt的频率最小值
		STRANDED_RATIO_THRES=0.1,     ## reads链偏性检测

	shell:
		"python {params.scriptPath}/step2_ntfreq_To_bigtables.py  -i  {input.configF} -r {input.refF}  -o {params.outputP}  \
		-D  {params.DEP_THRES} -A {params.MIN_DEP_THRES} -N {params.VALID_SIZE_THRES} -F {params.FREQ_THRES} -S {params.STRANDED_RATIO_THRES} "   #



rule convert_iSNVpytable_To_vcf :
	input:
		iSNVtable="/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/all.iSNV_with_SNP.pyResults.txt",
		ntfreqTable = "/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/out.samps-ntfreq.bigtable.txt",
		RefGnmF = "/shared/liuhj/HP/process/assembly/person_respectGnm/{person}.respect.fasta"

	output:
		outputP=directory("/shared/liuhj/HP/process/NGS_map2PersonResptGnm/table/{person}/vcf")
	params:
		scriptPath="/shared/liuhj/HP/scripts/HP-microevolution/scripts/iSNV_calling",
		MinFreq	= 0.02,
		#snv_MaxFreq	= 0.95,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_iSNVtable_To_vcf.py   -i {input.iSNVtable}  \
		-n {input.ntfreqTable}  -r  {input.RefGnmF} -o  {output.outputP} -m {params.MinFreq}  "








