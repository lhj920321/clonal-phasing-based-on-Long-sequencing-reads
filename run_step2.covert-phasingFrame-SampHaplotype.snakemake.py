###################################
#######phasing in iSNV peak  based on the NGS and CCS reliable posi and CCS alignment
###################################


##sample params
CCS_Samps=config["CCS_Sample"]
Reference_samp=config["REF_Sample"]
##SNV params
iSNVminFreq=config["iSNVminFreq"]
iSNVmaxFreq=config["iSNVmaxFreq"]
##path params
SNV_path=config["SNV_path"]
refFasta_path=config["refFasta_path"]
bam_path=config["bam_path"]
software_Path=config["software_Path"]
OUT_path=config["OUT_path"]
##out params
minPerct_covered=config["minPerct_covered"]
minhapNumNeed_forOut=config["minhapNumNeed_forOut"]
minDepth_forHapOut=config["minDepth_forHapOut"]
MaxdiffBaseNum_to_MergeReadsHaps=config["MaxdiffBaseNum_to_MergeReadsHaps"]




rule all:
	input:
	#1. convert frame
		expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}_2_{Ref_samp}.step1_convertColor.log",\
			samp={CCS_Samps},Ref_samp=Reference_samp,\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
		# expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}.Posi_AlleBases.txt",\
		# 	samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),

	#2. merge sample's haplotype
		expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}_2_{Ref_samp}.merge_phasing-Sample-Haps-mergePlot-HapStat.log",\
			samp={CCS_Samps},Ref_samp=Reference_samp,\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
		# expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/merge_phasing-Sample-Haps-mergePlot-HapStat/{samp}_2_{Ref_samp}.mergedHaplotypes.supportNum.txt",\
		# 	samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),



#Step step1 convert the phasing frame to the color block
rule revertDataFrame:
	input:
		#SampInfoF=expand("{software_Path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}.config.txt",\
		#	software_Path=software_Path,samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
		filt_posi_vcfF = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		refGnmP= config["refFasta_path"],
	output:
		RevertDataLog=temp(expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}_2_{Ref_samp}.step1_convertColor.log",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path)),
		Cites_AlleBasesF=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}.Posi_AlleBases.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
	params:
		scriptPath = config["software_Path"],
		RevertDataFrame_out_path = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		merge_out_path = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		Person=expand("{samp}_2_{Ref_samp}",samp={CCS_Samps},Ref_samp=Reference_samp),
		#samps=expand("{samp}",samp={CCS_Samps}),

	shell:
		"python  {params.scriptPath}/Step-step1-revert-dataFrame-Sample.py \
		-v  {input.filt_posi_vcfF}  -r {input.refGnmP}  -o  {params.RevertDataFrame_out_path} \
		-m {params.merge_out_path} -p {params.Person}; echo 'Step6' >>{output.RevertDataLog} "




#Step step2  merge sample block ,plot
rule Merge_and_Filt_Regions_to_plot_merge_filt:
	input:
		#SampInfoF=expand("{software_Path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}.config.txt",\
		#	software_Path=software_Path,samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
		filt_posi_vcfF = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		refGnmP=config["refFasta_path"],
		out_logF = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}_2_{Ref_samp}.step1_convertColor.log",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		Cites_AlleBasesF=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}.Posi_AlleBases.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
	output:
		RevertDataLog=temp(expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/{samp}_2_{Ref_samp}.merge_phasing-Sample-Haps-mergePlot-HapStat.log",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path)),
		outfile1 = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/merge_phasing-Sample-Haps-mergePlot-HapStat/{samp}_2_{Ref_samp}.mergedHaplotypes.supportNum.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
	params:
		scriptPath = config["software_Path"],
		RevertDataFrame_out_path = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/merge_phasing-Sample-Haps-mergePlot-HapStat",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		merge_out_path = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}",\
			samp={CCS_Samps},Ref_samp=Reference_samp,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,OUT_path=OUT_path),
		Person=expand("{samp}_2_{Ref_samp}",samp={CCS_Samps},Ref_samp=Reference_samp),
		minHapSuptNum = config["minDepth_forHapOut"],
		maxNumAlloToMerge = config["MaxdiffBaseNum_to_MergeReadsHaps"],
	shell:
		"python  {params.scriptPath}/Step-Step2-merge-SampleHaplotype.py  \
		-v {input.filt_posi_vcfF}  -r  {input.refGnmP} \
		-o  {params.RevertDataFrame_out_path} \
		-m {params.merge_out_path} -p {params.Person}  -S {params.minHapSuptNum}  \
		-P {params.maxNumAlloToMerge} -b {input.Cites_AlleBasesF}; echo 'Step7 ' >>{output.RevertDataLog} "


	






