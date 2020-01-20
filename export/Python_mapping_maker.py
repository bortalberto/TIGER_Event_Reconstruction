import sys
import numpy as np
import ROOT
import pickle

OS = sys.platform
if OS == 'win32':
    sep = '\\'
elif OS == 'linux':
    sep = '/'
else:
    print("ERROR: OS {} non compatible".format(OS))
    sys.exit()

class python_mapping_maker():
    def __init__(self):
        self.mapping_matrix = np.zeros((20,8,64), dtype=object)
    def make_map(self):
        """
        Extract the map in pytohn (map format for strip_display)
        :return:
        """
        # Create root file and tree
        map_file = ROOT.TFile.Open("mapping_IHEP_L2_2planari_penta.root")
        for event in map_file.tree:
            if event.pos_x>0:
                self.mapping_matrix[event.gemroc_id][event.SW_FEB_id][event.channel_id] = "X-{}".format(event.pos_x)
            elif event.pos_v>0:
                self.mapping_matrix[event.gemroc_id][event.SW_FEB_id][event.channel_id] = "V-{}".format(event.pos_v)
            else:
                self.mapping_matrix[event.gemroc_id][event.SW_FEB_id][event.channel_id] = "NaS"

        # print (self.baseline_matrix[6])

    def save_map(self):
        """
        Save the map in pickle encoded for Python 2
        :return:
            """
        with open("mapping.pickle", "wb+") as fileout:
                pickle.dump(self.mapping_matrix, fileout, protocol=2)
if __name__ == "__main__":
    mappatore=python_mapping_maker()
    mappatore.make_map()
    mappatore.save_map()
