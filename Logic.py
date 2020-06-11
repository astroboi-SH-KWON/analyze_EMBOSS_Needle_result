from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

import Util

class Logics:
    def __init__(self):
        self.tmp = ""

    def get_del_idx_seq(self, idx_key, del_dict, del_seq):
        if idx_key in del_dict:
            return self.get_del_idx_seq(idx_key + 1, del_dict, del_seq + del_dict[idx_key])
        else:
            return idx_key - 1, del_seq

    def get_sub_idx_seq(self, idx_key, sub_dict, seq_arr):
        sub_seq_from = seq_arr[0]
        sub_seq_to = seq_arr[1]
        if idx_key in sub_dict:
            sub_seq_from_char = sub_dict[idx_key].split("->")[0]
            sub_seq_to_char = sub_dict[idx_key].split("->")[1]
            return self.get_sub_idx_seq(idx_key + 1, sub_dict, [sub_seq_from + sub_seq_from_char, sub_seq_to + sub_seq_to_char])
        else:
            return idx_key - 1, sub_seq_from, sub_seq_to

    """
    :param
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
    """
    def get_sub_ins_del_list_dict_by_fnm(self, needle_dict):
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
        return result_dict

    def get_pairwise2_needle_dict(self, sources):
        result_dict = {}
        for i in range(len(sources)):
            tmp_list = []
            with open(sources[i], "r") as f:
                print(sources[i])
                print(f.readline())
                while True:
                    tmp_line = f.readline().replace("\n", "")
                    if tmp_line == "":
                        break
                    tmp_arr = tmp_line.split("\t")
                    idx = tmp_arr[0]
                    final_idx = tmp_arr[1]
                    ngs_read = tmp_arr[3]
                    ref_seq = tmp_arr[4]
                    ngs_read_needle, needle_result, ref_seq_needle = self.get_pairwise2_needle_result(ngs_read, ref_seq)
                    tmp_list.append([idx, final_idx, ngs_read_needle, needle_result, ref_seq_needle])

            result_dict[sources[i]] = tmp_list

        return result_dict

    """
    by using the BLOSUM62 matrix, together with a gap open penalty of 10 and a gap extension penalty of 0.5 (using globalds)
    """
    def get_pairwise2_needle_result(self, asequence, bsequence, matrx=blosum62, gap_open_penalty=10, extension_penalty=0.5):
        alignments = pairwise2.align.globalds(asequence.upper(), bsequence.upper(), matrx, -gap_open_penalty,
                                              -extension_penalty)
        align_arr = pairwise2.format_alignment(*alignments[0]).split("\n")
        return align_arr[0], align_arr[1], align_arr[2]



