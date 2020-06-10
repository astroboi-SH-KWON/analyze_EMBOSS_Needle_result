from Bio import SeqIO

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


