#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time,re 
from optparse import OptionParser
import pysam 




def readConfigF(config_F):
    config_ls = open(config_F,'r').readlines()
    for config_l in config_ls:
        if "ONT_PosiIndex_step" in config_l:
            ONT_PosiIndex_step = int(config_l.split("\n")[0].split("=")[1])
        #if "ONT_BGI_freq_maxDvalue" in config_l:
            #ONT_BGI_freq_maxDvalue = int(config_l.split("\n")[0].split("=")[1])

    return ONT_PosiIndex_step


def readRefGnm(ref_F):
    ref_seq = ''
    ref_F_ls = open(ref_F,'r').readlines()
    for ref_F_l in ref_F_ls:
        if '>'  not in ref_F_l:
            ref_seq = ref_seq + ref_F_l.split("\n")[0]

    return ref_seq



def readVcf(VcfF):
    CiteLst = {}
    for VcfFl in open(VcfF).readlines():
        if "#" not in VcfFl and VcfFl != "\n":
            Tags = VcfFl.split("\n")[0].split("\t")
            Cite = Tags[0]
            #Freq = float(Tags[9].split(":")[5].split("%")[0])
            CiteLst[int(Cite)] = ''

    return CiteLst               



def ONTreads_SmallAltposi_char(POS,SEQ,CIGAR,QNAME,SEQ_len,ref_seq,bam_reads_line,posi_list,out_smallAltPosi_F_O,out_line_tongji_dict):
    reads_base_count = 0
    new_seq_list = []
    new_seq = ''

    SEQ_list = list(SEQ)      
    for number_index in range(0,len(CIGAR)):        
        CIGAR_char = CIGAR[number_index][0]
        CIGAR_char_num = CIGAR[number_index][1]
        left = reads_base_count
        if CIGAR_char == 4 or CIGAR_char == 0 or CIGAR_char == 1 :   ###4--S,0--M,1--I
            reads_base_count = reads_base_count + int(CIGAR_char_num)     ##当有S、M、I出现时,因为这些是在reads确实存在的碱基。因此，         
        if CIGAR_char == 0:     ##将M的所有字符从reads上截取下来，放到新的seq中
            new_seq = new_seq + SEQ[left:reads_base_count]
        elif CIGAR_char == 2:     ##将D 对应的数字，表示成相应个数的‘-’  其中2 --- D
            for kong_num in range(0,int(CIGAR_char_num)):
                new_seq = new_seq + "-"
                
    new_seq_len = len(new_seq)    
    read_alt_list = []
    out_line_1 = ''
    out_line = str(QNAME) + "\t" + str(SEQ_len) + "\t"
    #out_line_Lst =
    for aim_posi in posi_list:
        if aim_posi >= POS  and aim_posi <  POS + new_seq_len:
            aim_posi_read_char = new_seq [aim_posi - POS - 1 ]
            aim_posi_ref_char = ref_seq [ aim_posi - POS - 1 ]
            
            #if (aim_posi_read_char != '-') and (aim_posi_read_char not in  BGIaimPosi_altFreq_dic[str(aim_posi)].keys()) and (aim_posi_read_char not in  BGIaimPosi_refFreq_dic[str(aim_posi)].keys()) :
                #aim_posi_read_char = aim_posi_read_char + '*'
            out_line =  out_line    + str(aim_posi)  + aim_posi_read_char + "_"
            out_line_1 = out_line_1  + str(aim_posi)  + aim_posi_read_char + "_"


    out_smallAltPosi_F_O.write(out_line + "\n")

    if out_line_1 != '':
        if out_line_1 not in out_line_tongji_dict.keys():
            out_line_tongji_dict[out_line_1] = 0
            out_line_tongji_dict[out_line_1] = out_line_tongji_dict[out_line_1]  + 1
        else:
            out_line_tongji_dict[out_line_1] = out_line_tongji_dict[out_line_1]  + 1

    



def main():

    ################################################################################################################################
    # program
    ################################################################################################################################
    startTime = time.time()
    #ONT_PosiIndex_step = readConfigF(config_F)

    small_altPosi = readVcf(snv_vcf_F)
    small_altPosi_lst = list(small_altPosi)
    small_altPosi_lst.sort()

    ref_seq = readRefGnm(ref_F)    

    #####------------------------------------------------------------------------------
    ###########  将筛选后剩下的可靠位点，建立posi的索引
    #ONT_PosiIndex_step = 20000    ###参数 ：建立位点索引的范围大小
    for left_limit in range(0,len(ref_seq),ONT_PosiIndex_step):
        #print left_limit
        right_limit = left_limit + ONT_PosiIndex_step
        locals()['limit_list_%s'%left_limit] = []
        for aim_posi in small_altPosi_lst:  
            if aim_posi >= left_limit  and aim_posi < right_limit:
                locals()['limit_list_%s'%left_limit].append(aim_posi)
          
    out_line_tongji_dict = {}


    #####重新读取ONT bam文件，提取出筛选得到可靠位点以后，这些可靠位点在每条reads 上对应的组成
    out_smallAltPosi_F_O = open(FiltPosi_ONTsnvsType_F,'a')
    bf = pysam.AlignmentFile(ONT_bam_F, 'rb')
    bam_line_count = 0
    for bam_line in bf:
        bam_line_count = bam_line_count + 1
        if bam_line_count % 10000 == 0:
            print("reliably posi        bam_line:" + str(bam_line_count))
        POS = bam_line.pos
        CIGAR = bam_line.cigar
        SEQ = bam_line.seq
        QNAME = bam_line.qname
        SEQ_len = len(SEQ)
        #print(str(POS) + "\t" + str(len(SEQ)))
        left_flag = int(POS/ONT_PosiIndex_step)*ONT_PosiIndex_step
        right_Posi = POS + len(SEQ)
        if right_Posi  > len(ref_seq):
            right_Posi = len(ref_seq)
        if (POS >= left_flag)  and  (right_Posi < left_flag + ONT_PosiIndex_step ):
            aim_posi_list = locals()['limit_list_%s'%left_flag]                 
        else:
            aim_posi_list = [] 
            if POS > 0:               
                for FLAG in range(int(left_flag), (int(right_Posi/ONT_PosiIndex_step) + 1)*ONT_PosiIndex_step,ONT_PosiIndex_step):   #int(fanwei_flag) +    
                    aim_posi_list = aim_posi_list + locals()['limit_list_%s'%FLAG]    # + locals()['limit_list_%s'%fanwei_flag_next]
        ONTreads_SmallAltposi_char(POS,SEQ,CIGAR,QNAME,SEQ_len,ref_seq,bam_line,aim_posi_list,out_smallAltPosi_F_O,out_line_tongji_dict)


    print("PASS :" + str(len(small_altPosi_lst)) )



    ##sort 
    #for BGI_aim in BGI_snv_posi_list
    FiltPosi_ONTsnvsType_F
    out_ONT_tongji_F_O = open(FiltPosi_ONTsnvsType_stat_F,'a')

    for out_line_tongji in out_line_tongji_dict.keys():
            #if str("_" + str(BGI_aim)) in out_line_tongji:  
        genotype_lst = out_line_tongji.split("\t")
        first_posi = re.findall('\d+',genotype_lst[0])[0]
        last_posi = re.findall('\d+',genotype_lst[0])[-1]
    

      
        out_line = str(out_line_tongji_dict[out_line_tongji]) + "\t" + "S-" + first_posi + "-" + last_posi + "\t" + out_line_tongji
        out_ONT_tongji_F_O.write(out_line + "\n")            
                #del out_line_tongji_dict[out_line_tongji]

    #out_ONTreads_snvsType_FO.close()
    #out_smallAltPosi_F_O.close()
    #out_freq_F_O.close()
    endTime = time.time()
    sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")




if __name__ == "__main__":
################################################################################################################################
# Parameters
################################################################################################################################
    usage = "usage:python  %prog [option]"
    parser = OptionParser(usage=usage)
    #parser.add_option("-c","--config_F",
                        #dest = "config_F",
                        #default = "",
                        #metavar = "file",
                        #help = "Config file [required]")
    parser.add_option("-r","--ref_F",
                        dest = "ref_F",
                        default = "",
                        metavar = "string",
                        help = "reference genome file [required]")
    parser.add_option("-b","--ONT_bam_F",
                        dest = "ONT_bam_F",
                        default = "",
                        metavar = "file",
                        help = "bam file  [required]")
    parser.add_option("-v","--snv_vcf_F",
                        dest = "snv_vcf_F",
                        default = "",
                        metavar = "file",
                        help = "snv -- vcf fomated file. [required]")
    parser.add_option("-o","--out_Path",
                        dest = "out_Path",
                        default = "",
                        metavar = "path",
                        help = "Output Path. [required]")
    parser.add_option("-s","--ONT_PosiIndex_step",
                        dest = "ONT_PosiIndex_step",
                        default = "200",
                        metavar = "int",
                        help = "step to build the ONT Posi Index  . [required]")

    
    (options,args) = parser.parse_args()
    #config_F       = os.path.abspath(options.config_F)
    ref_F          = os.path.abspath(options.ref_F)
    ONT_bam_F      = os.path.abspath(options.ONT_bam_F)
    snv_vcf_F      = os.path.abspath(options.snv_vcf_F)
    out_Path       = os.path.abspath(options.out_Path)
    ONT_PosiIndex_step = int(options.ONT_PosiIndex_step)
    #out_ONTreads_snvsType_F = out_Path + "/1.out-bam_SnvType.txt"
    #out_filtePosi_F         = out_Path + "/1.out-Filted_Cites.txt"
    FiltPosi_ONTsnvsType_F  = out_Path + "/1.out-SNVtype_of_Filted_Cites.txt"
    print(FiltPosi_ONTsnvsType_F)
    FiltPosi_ONTsnvsType_stat_F = out_Path + "/1.out-SNVtype_of_Filted_Cites.Stat.txt"
    print(FiltPosi_ONTsnvsType_stat_F)
    #FiltPosi_lines_F = out_Path + "/1.out-Filted_reliablePosi.txt"


    ################################################################################################################################
    # Check input files
    ################################################################################################################################
    if (not os.path.exists(out_Path)):
        os.mkdir(out_Path)
    #if (not os.path.exists(config_F)):
        #print("config")

    if (not os.path.exists(ONT_bam_F)):
        print("long reads bam file exists !")
    if (not os.path.exists(snv_vcf_F)):
        print("reliable iSNV vcf format file exists !")
    #if (os.path.exists(out_ONTreads_snvsType_F)) :
        #os.remove(out_ONTreads_snvsType_F)
    #if (os.path.exists(out_filtePosi_F)) :
        #os.remove(out_filtePosi_F)
    #if (os.path.exists(FiltPosi_ONTsnvsType_F)) :
        #os.remove(FiltPosi_ONTsnvsType_F)
    #if (os.path.exists(FiltPosi_ONTsnvsType_stat_F)) :
        #os.remove(FiltPosi_ONTsnvsType_stat_F)


    main()




