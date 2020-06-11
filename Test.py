from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62

import pandas as pd
import Bio as bio
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
INPUT_DIR = "input/*test.txt"
############### end setting env ################

def main():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + SOURCE_DIR)

    needle_dict = util.get_result_dict_by_fnm(sources)

    result_dict = logic.get_sub_ins_del_list_dict_by_fnm(needle_dict)

    util.make_excel(WORK_DIR + "excel_output/result_", result_dict)

def pairwise2_main():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + INPUT_DIR)

    needle_dict = logic.get_pairwise2_needle_dict(sources)

    result_dict = logic.get_sub_ins_del_list_dict_by_fnm(needle_dict)

    util.make_excel(WORK_DIR + "excel_output/pairwise2_result_", result_dict)

def pairwise2_test():
    asequence = "TATATATCTTGTGGAAAGGACGAAACACCGGAAGCTGTACTTCAAAAAAGTTAGTACATTTTTTTCATATCTGCACTCACTCTCTGCTGAAGCTGTACTTCAAAAAATGGATGACATGAAGAAGATAGCTTGGCGTACCGCGATCTCTACTCTACCCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA"
    bsequence = "TATATATCTTGTGGAAAGGACGAAACACCGGAAGCTGTACTTCAAAAAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCATCgATTTTTTGAAGTACATTTTTTTCATATCTGCACTCACTCTCTGCTGAAGCTGTACTTCAAAAAATGGATGACATGAAGAAGATAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA"
    alignments = pairwise2.align.globalxx(asequence, bsequence)
    # print(alignments[0][0])
    # print(alignments[0][1])
    # print("")

    for tmp_tup in alignments:
        # print(tmp_tup[0])
        print(tmp_tup[1])
        # print("")

def pairwise2_w_opt_test():
    asequence = "TATATATCTTGTGGAAAGGACGAAACACCGGAAGCTGTACTTCAAAAAAGTTAGTACATTTTTTTCATATCTGCACTCACTCTCTGCTGAAGCTGTACTTCAAAAAATGGATGACATGAAGAAGATAGCTTGGCGTACCGCGATCTCTACTCTACCCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA"
    bsequence = "TATATATCTTGTGGAAAGGACGAAACACCGGAAGCTGTACTTCAAAAAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCATCgATTTTTTGAAGTACATTTTTTTCATATCTGCACTCACTCTCTGCTGAAGCTGTACTTCAAAAAATGGATGACATGAAGAAGATAGCTTGGCGTAcCgcGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA"
    # by using the BLOSUM62
    # matrix, together with a gap open penalty of 10 and a gap extension penalty of 0.5 (using globalds)
    gap_open_penalty = 10
    extension_penalty = 0.5
    alignments = pairwise2.align.globalds(asequence.upper(), bsequence.upper(), blosum62, -gap_open_penalty, -extension_penalty)
    # for tmp_tup in alignments:
    #     print(tmp_tup[0])
    #     print(tmp_tup[1])
    #     print("")

    print(pairwise2.format_alignment(*alignments[0]).split("\n"))
    print(len(pairwise2.format_alignment(*alignments[0]).split("\n")[0]))
    print(len(pairwise2.format_alignment(*alignments[0]).split("\n")[1]))
    print(len(pairwise2.format_alignment(*alignments[0]).split("\n")[2]))




start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# main()
# pairwise2_w_opt_test()
pairwise2_main()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))