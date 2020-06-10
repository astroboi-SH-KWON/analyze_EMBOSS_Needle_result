import pandas as pd
import Bio as bio
from Bio import SeqIO
from time import clock
import glob
import multiprocessing as mp
from threading import Thread
import numpy as np


import Util
import Logic
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/YuGooSang_KimHuiKwon/20200609/WORK_DIR/"

SOURCE_DIR = "crawler_output/*test.txt"
############### end setting env ################

def main():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + SOURCE_DIR)

    needle_dict = util.get_result_dict_by_fnm(sources)

    result_dict = logic.get_sub_ins_del_list_dict_by_fnm(needle_dict)

    util.make_excel(WORK_DIR + "excel_output/result_", result_dict)

def test():
    util = Util.Utils()

    needle_dict = {'D:/000_WORK/YuGooSang_KimHuiKwon/20200609/WORK_DIR/first_excel_output\\result_gDNA_0609.txt':
        [['1', 'Group1,2_RT/20-PBS/7-#Target723'
        , 'AATATATCTTGTGGAAAGGACGAAACACCG--CATACTCGGGCGC-------CGGGGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGACCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACATGCCAGGTGGACGAGTTTTCTTGCTTTTTTTGATACTCTGTCTGTACTACAACGCCCATTTCCGCAAGAAAACTGGTCTACCTGGCATGTTCAGCTTGGCGTACCGCGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA'
        , '.|||||||||||||||||||||||||||||  ||| .||   |||       |     ||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
        , 'TATATATCTTGTGGAAAGGACGAAACACCGCCCAT-TTC---CGCAAGAAAAC-----GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACATGCCAGGTAGACgAGTTTTCTTGCTTTTTTTGATACTCTGTCTGTACTACAACGCCCATTTCCGCAAGAAAACTGGTCTACCTGGCATGTTCAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA', '293', '293', '293', 'O']
        , ['2', 'Group1,2_RT/12-PBS/11-#Target1948'
        , 'TATATATCTTGTGGAAAGGACGAAACACCGAAGTCCGTCAGATTCTATCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCATACCACGAGATAGAATCTGACGTTTTTTTCGTACTCATATATACATATCTCTAAGTCCGTCAGATTCTATCTGGTGGTATCTCCAGGTGAAGCTTGGCGTACCGCGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA'
        , '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
        , 'TATATATCTTGTGGAAAGGACGAAACACCGAAGTCCGTCAGATTCTATCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCATACCACgAGATAGAATCTGACGTTTTTTTCGTACTCATATATACATATCTCTAAGTCCGTCAGATTCTATCTGGTGGTATCTCCAGGTGAAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA',  '280', '280', '280', 'O']
        , ['3', 'Group1,2_RT/20-PBS/17-#Target833'
        , '-----------------------------------------------------------------------------------------CG-------CTTGAAAAAGTGGCACCGAGTCGGTGCTTACCTCTTTGGATCGTGATCACAATCCTCCAGATGCTTTTTTTCAGATAGCATACTGTATACTGGGCATCTGGAGGATTGTGATCAGGATCCAAAGAGGTAATGAGCTTGGCGTACCGCGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAACGTGCACGTGACACGTTCCAGACCGTACATGCTTACATGGGATGAAGCTTGGCGTAACTAGATCTTGAGACAAATGGCAGTATT'
        , '                                                                                         ||       ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                                                                    '
        , 'TATATATCTTGTGGAAAGGACGAAACACCGCATCTGGAGGATTGTGATCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTACCTCTTTGGATCgTGATCACAATCCTCCAGATGCTTTTTTTCAGATAGCATACTGTATACTGGGCATCTGGAGGATTGTGATCAGGATCCAAAGAGGTAATGAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA------------------------------------------------------------------------------------', '378', '378', '378', 'O']
        , ['4', 'Group1,2_RT/12-PBS/9-#Target489'
        , 'TATATATCTTGTGGAAAGGACGAAACACCGGCGCGGAACAGGTCG--ATC-TGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAAGTACCGTTTGATGCCGCTGTTTTTTTCATACACGACACACATCTGAGGTCGTTCACCAGCGGCATCAAAGGGTACTTCATGGCGCATAGCTTGGTGTACCGCGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA'
        , '||||||||||||||||||||||||||||| |||...|.|||  ||  ||| .||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
        , 'TATATATCTTGTGGAAAGGACGAAACACC-GCGTTCACCAG--CGGCATCAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAAGTACCgTTTGATGCCGCTGTTTTTTTCATACACGACACACATCTGAGGTCGTTCACCAGCGGCATCAAAGGGTACTTCATGGCGCATAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA', '281', '281', '281', 'O']
        ]}

    result_dict = {}
    for fnm_key, val_list in needle_dict.items():
        result_dict.update({fnm_key: []})
        for val_arr in val_list:
            final_index = val_arr[1]
            ngs_read = val_arr[2].upper()
            ref_seq = val_arr[4].upper()
            re_idx = 0
            sub_dict = {}
            ins_dict = {}
            del_dict = {}
            for i in range(len(ngs_read)):
                if ngs_read[i] != ref_seq[i]:
                    if ngs_read[i] == "-":
                        re_idx += 1
                        del_dict.update({re_idx: ref_seq[i]})
                    elif ref_seq[i] == "-":
                        if re_idx in ins_dict:
                            ins_dict[re_idx] += ngs_read[i]
                        else:
                            ins_dict.update({re_idx: ngs_read[i]})
                    else:
                        re_idx += 1
                        sub_dict.update({re_idx: ref_seq[i] + "->" + ngs_read[i]})
                else:
                    re_idx += 1
            result_dict[fnm_key].append([final_index, sub_dict, ins_dict, del_dict, re_idx])

    util.make_excel(WORK_DIR + "excel_output/result_", result_dict)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
main()
# test()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))