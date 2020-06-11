from time import clock

import Util
import Logic
############### start to set env ################
WORK_DIR = "D:/000_WORK/YuGooSang_KimHuiKwon/20200609/WORK_DIR/"

SOURCE_DIR = "crawler_output/*.txt"
INPUT_DIR = "input/*.txt"
############### end setting env ################

def main():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + SOURCE_DIR)

    needle_dict = util.get_result_dict_by_fnm(sources)

    result_dict = logic.get_sub_ins_del_list_dict_by_fnm(needle_dict)

    util.make_excel(WORK_DIR + "excel_output/result_time", result_dict)

def pairwise2_main():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + INPUT_DIR)

    needle_dict = logic.get_pairwise2_needle_dict(sources)

    result_dict = logic.get_sub_ins_del_list_dict_by_fnm(needle_dict)

    util.make_excel(WORK_DIR + "excel_output/pairwise2_result_", result_dict)

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# main()
pairwise2_main()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))