import data as dt
import glob2
import glob
import sys
import numpy as np
import os
import ROOT
from ROOT import gROOT, AddressOf
import matplotlib.pyplot as plt
import pickle
OS = sys.platform
if OS == 'win32':
    sep = '\\'
elif OS == 'linux2':
    sep = '/'
else:
    print("ERROR: OS {} non compatible".format(OS))
    sys.exit()


class baseline_importer():
    
    def __init__(self, scan_path):
        self.run_path = scan_path
        self.baseline = {}
        self.baseline_from_calib_matrix=np.zeros((20, 8, 64))
        self.real_baseline_T = np.zeros((20, 8, 64))
        self.real_baseline_E = np.zeros((20, 8, 64))
        self.thr_matrix_T = np.zeros((20, 8, 64))
        self.thr_matrix_E = np.zeros((20, 8, 64))
        self.distance_from_baseline = np.zeros((20,8,64))
        self.diff_matrix = np.zeros((20,8,3))
        self.primo_giro=True
        
    def extract_baseline(self):
        """
        Excrat the baseline values and write them in a pickle (and in a root file)
        :return:
        """
        # Create root file and tree
        map_file = ROOT.TFile.Open("mapping_IHEP_L2_2planari_penta.root")
        mapping_matrix = {1: {}, 2: {}, 3: {},0:{}}
        for event in map_file.tree:
            mapping_matrix[event.layer_id][event.HW_FEB_id, event.chip_id] = [event.gemroc_id,event.SW_FEB_id]
        gROOT.ProcessLine('struct TreeStruct {\
                int layer_id;\
                int hardware_FEB_id;\
                int chip_id;\
                int channel_id;\
                float baseline;\
                };')
        rootFile = ROOT.TFile("baselines.root", 'recreate')
        tree = ROOT.TTree('tree', '')
        mystruct = ROOT.TreeStruct()
        for key in ROOT.TreeStruct.__dict__.keys():
            if '__' not in key:
                formstring = '/F'
                if isinstance(mystruct.__getattribute__(key), int):
                    formstring = '/I'
                tree.Branch(key, AddressOf(mystruct, key), key + formstring)
        for root, dirs, files in os.walk(self.run_path ):
            for scan_file, (layer, HW_FEB, chip, others) in glob2.iglob(root+ sep + "L*FEB*_c*_VTH*.dat", with_matches=True,recursive=True):
                importer = dt.Import_Data("./", scan_file, False)
                try:
                    GEMROC_ID, TIGER_ID = mapping_matrix[int(layer)][int(HW_FEB), int (chip)]
                    for ch in range(0, 64):
                        data = importer.labVIEW_data[ch][0]['tot_evts']
                        self.baseline_from_calib_matrix[GEMROC_ID][TIGER_ID][ch]=int((np.argmax(data)))
                        # if GEMROC_ID==6 and TIGER_ID==7:
                        #     input("Ho preso per gemroc {}, tiger {} dal file {}".format(GEMROC_ID,TIGER_ID,scan_file))
                        mystruct.layer_id = int(layer)
                        mystruct.hardware_FEB_id = int(HW_FEB)
                        mystruct.chip_id = int(chip)
                        mystruct.channel_id = int(ch)
                        mystruct.baseline = int((np.argmax(data)))
                        tree.Fill()
                except KeyError:
                    print ("{} is spare FEB".format(HW_FEB))
        rootFile.Write()
        rootFile.Close()
        # print (self.baseline_matrix[6])

    def extract_baseline_L1(self):
        """
        Excrat the baseline values and write them in a pickle (and in a root file)
        :return:
        """
        # Create root file and tree
        map_file = ROOT.TFile.Open("mapping_IHEP_L2_2planari_penta.root")
        mapping_matrix = {1: {}, 2: {}, 3: {},0:{}}
        for event in map_file.tree:
            mapping_matrix[event.layer_id][event.HW_FEB_id, event.chip_id] = [event.gemroc_id,event.SW_FEB_id]
        print( mapping_matrix[1][10,2])
        T_list=[]
        E_list=[]
        T_mapped={}
        E_mapped={}


        with open("L1/Baseline_T_1.pickle", "rb") as file_sistema:
            T_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_E_1.pickle", "rb") as file_sistema:
            E_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_T_2.pickle", "rb") as file_sistema:
            T_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_E_2.pickle", "rb") as file_sistema:
            E_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_T_3.pickle", "rb") as file_sistema:
            T_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_E_3.pickle", "rb") as file_sistema:
            E_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_T_4.pickle", "rb") as file_sistema:
            T_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_E_4.pickle", "rb") as file_sistema:
            E_list.append(pickle.load(file_sistema,encoding="latin1"))
        with open("L1/Baseline_T_5.pickle", "rb") as file_sistema:
            T_list.append(pickle.load(file_sistema, encoding="latin1"))
        with open("L1/Baseline_E_5.pickle", "rb") as file_sistema:
            E_list.append(pickle.load(file_sistema, encoding="latin1"))

        RUN_map=[]
        RUN_map.append({
            0: (5, 1),
            1: (5, 2),
            2: (2, 1),
            3: (2, 2),
            4: (13, 1),
            5: (13, 2),
            6: (10, 1),
            7: (10, 2)

        })
        RUN_map.append({
            0: (24, 1),
            1: (24, 2),
            2: (25, 1),
            3: (25, 2),
            4: (22, 1),
            5: (22, 2),
            6: (21, 1),
            7: (21, 2)

        })
        RUN_map.append({
            0: (4, 1),
            1: (4, 2),
            2: (26, 1),
            3: (26, 2),
            4: (19, 1),
            5: (19, 2),
            6: (29, 1),
            7: (29, 2)

        })
        RUN_map.append({
            0: (30, 1),
            1: (30, 2),
            2: (34, 1),
            3: (34, 2),
            4: (32, 1),
            5: (32, 2),
            6: (61, 1),
            7: (61, 2)

        })
        RUN_map.append({
            0: (8, 1),
            1: (8, 2),
            2: (33, 1),
            3: (33, 2),
            4: (27, 1),
            5: (27, 2),
            6: (16, 1),
            7: (16, 2)

        })

        for run_map,baseline_t,baseline_e in zip(RUN_map,T_list,E_list): ##Builds the mapped version of the baseline dicts
            for key in run_map:
                T_mapped[run_map[key][0],run_map[key][1]]=baseline_t["GEMROC 0"]["TIG{}".format(key)]
                E_mapped[run_map[key][0],run_map[key][1]]=baseline_e["GEMROC 0"]["TIG{}".format(key)]

        for key in mapping_matrix[1]:
            for channel in range (0,64):
                self.real_baseline_T[mapping_matrix[1][key][0]][mapping_matrix[1][key][1]][channel] = T_mapped[key[0],key[1]]["CH{}".format(channel)][3]-5
                self.real_baseline_E[mapping_matrix[1][key][0]][mapping_matrix[1][key][1]][channel] = E_mapped[key[0],key[1]]["CH{}".format(channel)][3]-4

    def save_in_pickle(self):
        with open("baseline_from_calib","wb+") as fileout:
            pickle.dump(self.baseline_from_calib_matrix, fileout)

    def save_in_pickle_real_baseline(self):
        """
        Save the real baselines in a pickle file
        :return:
        """
        with open("real_baseline_T.pickle", "wb+") as fileout:
            pickle.dump(self.real_baseline_T, fileout, protocol=2)  # Protocol=2 for python 2

        with open("real_baseline_E.pickle", "wb+") as fileout:
            pickle.dump(self.real_baseline_E, fileout, protocol=2)  # Protocol=2 for python 2

    def load_real_baseline(self):
        with open("real_baseline_T.pickle", "rb") as filein:
            self.real_baseline_T = pickle.load(filein)  # Protocol=2 for python 2

        with open("real_baseline_E.pickle", "rb") as filein:
            self.real_baseline_E = pickle.load(filein)  # Protocol=2 for python 2

    def confront_baselines(self):
        """
        Confront baselines from calibration and actual baseline on the detector
        :return:
        """
        with open("baseline_from_calib","rb") as file_calib:
            bas_from_calib = pickle.load(file_calib)
        with open("Baseline_T.pickle", "rb") as file_sistema:
            self.bas_from_sys = pickle.load(file_sistema,encoding="latin1")
        with open("Baseline_E.pickle", "rb") as file_sistema:
            self.bas_from_sys_E = pickle.load(file_sistema, encoding="latin1")
        for GEMROC in range (4,11): #For layer 2
            print ("GEMROC {}".format(GEMROC))
            for TIGER in range (0,8):
                print("TIGER {}".format(TIGER))
                diff_list = []
                for channel in range (61,64):
                    """ Use those channels as reference"""
                    if bas_from_calib[GEMROC][TIGER][channel]==0:
                        print ("Error, no scan baseline for {} {}".format(GEMROC,TIGER))
                    print(bas_from_calib[GEMROC][TIGER][channel])
                    # diff_list.append((bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]-bas_from_calib[GEMROC][TIGER][channel]))
                    # print ((bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]-bas_from_calib[GEMROC][TIGER][channel]))
                    diff_list.append(((self.bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]-bas_from_calib[GEMROC][TIGER][channel])))
                    self.diff_matrix[GEMROC][TIGER][channel-61] = self.bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]-bas_from_calib[GEMROC][TIGER][channel]
                average_diff = (np.average(diff_list))
                for channel in range(0,64):
                    # print("Original {}".format(bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]))
                    # print(round(bas_from_sys["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][2]+average_diff))
                    self.real_baseline_T[GEMROC][TIGER][channel] = (round(bas_from_calib[GEMROC][TIGER][channel] + average_diff))
                    self.real_baseline_E[GEMROC][TIGER][channel] = (round(bas_from_calib[GEMROC][TIGER][channel] + average_diff))

                #     try:
                #         self.real_baseline_E[GEMROC][TIGER][channel] = (round(self.bas_from_sys_E["GEMROC {}".format(GEMROC)]["TIG{}".format(TIGER)]["CH{}".format(channel)][1] - average_diff))
                #     except Exception as e:
                #         print (e)
                #     plt.plot(diff_list)
                #     plt.show()
    def import_run_thr(self,filename):
        """
        Funcion to import the channel configuration from old RUNS. Needs the pickle containing the data
        :return:
        """
        File_name = filename

        with open(File_name, 'rb') as f:
            old_conf_dict = pickle.load(f)

        for GEMROC_key, dict in old_conf_dict.items():
            GEMROC_id=GEMROC_key.split(" ")[1]
            for TIGER_key, dict2 in dict.items():
                if TIGER_key.split(" ")[0] == "TIGER":
                    TIGER_id = int(TIGER_key.split(" ")[1])
                    for channel_key, dict3 in dict2.items():
                        if channel_key.split(" ")[0] == "Ch":
                            channel_id = int(channel_key.split(" ")[1])
                            # print ("{} - {} - {}".format(GEMROC_id, TIGER_id, channel_id))
                            try:
                                self.thr_matrix_T[int(GEMROC_id)][int(TIGER_id)][int(channel_id)]=dict3["Vth_T1"]
                                self.thr_matrix_E[int(GEMROC_id)][int(TIGER_id)][int(channel_id)]=dict3["Vth_T2"]
                            except KeyError as e:
                                print ('I a KeyError - missing %s. Probably a GEMROC is offline '% str(e))
                                break
                            # print self.GEMROC_reading_dict[GEMROC_key].c_inst.Channel_cfg_list[TIGER_id][channel_id]

    def save_delta_VTHR_root(self):
        """
        Save the information about thr and baseline in a root file
        :return:
        """

        if self.primo_giro:
            gROOT.ProcessLine('struct TreeStruct2 {\
                            int layer_id;\
                            int gemroc_id;\
                            int software_feb_id;\
                            int channel_id;\
                            int baseline_t;\
                            int baseline_e;\
                            int vth1;\
                            int vth2;\
                            int delta_vth1_baseline;\
                            int delta_vth2_baseline;\
                            float thr_t_fC;\
                            float thr_e_fC;\
            };')
            self.primo_giro=False
        rootFile = ROOT.TFile("delta_vth.root", 'recreate')
        tree = ROOT.TTree('tree', '')
        mystruct = ROOT.TreeStruct2()
        for key in ROOT.TreeStruct2.__dict__.keys():
            if '__' not in key:
                formstring = '/F'
                if isinstance(mystruct.__getattribute__(key), int):
                    formstring = '/I'
                tree.Branch(key, AddressOf(mystruct, key), key + formstring)

        for GEMROC in range(0, 11):
            if GEMROC < 4:
                lyr = 1
            elif GEMROC < 11:
                lyr = 2
            else:
                lyr = 3
            for TIGER in range(0, 8):
                for ch in range(0, 64):
                    mystruct.layer_id = int(lyr)
                    mystruct.gemroc_id = (int(GEMROC))
                    mystruct.software_feb_id = int(TIGER)
                    mystruct.channel_id = int(ch)
                    mystruct.baseline_t = int(self.real_baseline_T[GEMROC][TIGER][ch])
                    mystruct.baseline_e = int(self.real_baseline_E[GEMROC][TIGER][ch])
                    mystruct.vth1 = int(self.thr_matrix_T[GEMROC][TIGER][ch])
                    mystruct.vth2 = int(self.thr_matrix_E[GEMROC][TIGER][ch])
                    mystruct.delta_vth1_baseline = int(self.real_baseline_T[GEMROC][TIGER][ch]-self.thr_matrix_T[GEMROC][TIGER][ch])
                    mystruct.delta_vth2_baseline = int(self.real_baseline_E[GEMROC][TIGER][ch]-self.thr_matrix_E[GEMROC][TIGER][ch])
                    mystruct.thr_t_fC = float(self.real_baseline_T[GEMROC][TIGER][ch]-self.thr_matrix_T[GEMROC][TIGER][ch])*0.42
                    mystruct.thr_e_fC = float(self.real_baseline_E[GEMROC][TIGER][ch]-self.thr_matrix_E[GEMROC][TIGER][ch])*0.42
                    tree.Fill()

        rootFile.Write()
        rootFile.Close()
        os.system("root -l convert_VTH.cxx")
        os.system("rm delta_vth.root")
if __name__ == "__main__":
    importatore = baseline_importer("./scans")
    for runfile,(run_numb,) in glob2.iglob("conf_files/CONF_run_*.pkl", with_matches=True):
        # importatore.extract_baseline_L1()
        # importatore.extract_baseline()
        # importatore.save_in_pickle()
        # importatore.confront_baselines()
        # importatore.save_in_pickle_real_baseline()
        importatore.import_run_thr(runfile)
        importatore.load_real_baseline()
        importatore.save_delta_VTHR_root()
        os.system("mv delta_vth_mapped.root thr_run_{}.root".format(run_numb))
