import binascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import ROOT
import sys
from ROOT import gROOT, AddressOf



class reader:
    def __init__(self, GEMROC_ID, MODE, RUN):
        self.MODE = int(MODE)
        self.GEMROC_ID = int(GEMROC_ID)
        self.RUN=int(RUN)
        self.thr_scan_matrix = np.zeros((8, 64))  # Tiger,Channel
        self.thr_scan_frames = np.ones(8)
        self.thr_scan_rate = np.zeros((8, 64))

    def __del__(self):
        print ("Done")

    def write_root(self, path):
        gROOT.ProcessLine('struct TreeStruct {\
                int runNo;\
                int layer_id;\
                int gemroc_id;\
                int tiger_id;\
                int channel_id;\
                int tac_id;\
                float tcoarse;\
                float tcoarse_10b;\
                float ecoarse;\
                float tfine;\
                float efine;\
                float charge_SH;\
                int count;\
                int count_ori;\
                int count_new;\
                int timestamp;\
                int l1ts_min_tcoarse;\
                int lasttigerframenum;\
                int local_l1_count;\
                int count_mismatch;\
                float delta_coarse;\
                };')
        from ROOT import TreeStruct

        rname = path.replace(".dat",".root")

        rootFile = ROOT.TFile(rname, 'recreate')
        tree = ROOT.TTree('tree', '')
        mystruct = TreeStruct()

        for key in TreeStruct.__dict__.keys():
            if '__' not in key:
                formstring = '/F'
                if isinstance(mystruct.__getattribute__(key), int):
                    formstring = '/I'
                tree.Branch(key, AddressOf(mystruct, key), key + formstring)


        statinfo = os.stat(path)
        packet_header = -1
        packet_tailer = -1
        l1count = -1

        l1count_new=[]
        lschannel_id=[]
        lstac_id=[]
        lstcoarse=[]
        lstcoarse_10b=[]
        lsecoarse=[]
        lstfine=[]
        lsefine=[]
        lscharge_SH=[]
        lstigerid=[]
        lsl1ts_min_tcoarse=[]
        lslasttigerframenum=[]
        lscount_mismatch=[]
        lsdelta_coarse=[]
        l1timestamp = -1
        gemrocid = -1

        pre_timestamp = 0

        hitcounter=0
        max_hitcount=1000000000000
        flag_swap1=False
        flag_swap2=False

        with open(path, 'rb') as f:
            for i in range(0, statinfo.st_size // 8):
                    data = f.read(8)
                    #hexdata = binascii.hexlify(data)
                    if sys.version_info[0]== 2:
                        hexdata = str(binascii.hexlify(data))
                    else:
                        hexdata = str(binascii.hexlify(data), 'ascii')
                    string= "{:064b}".format(int(hexdata,16))
                    inverted=[]
                    for i in range (8,0,-1):
                        inverted.append(string[(i-1)*8:i*8])
                    string_inv="".join(inverted)
                    int_x = int(string_inv,2)

                # for x in range(0, len(hexdata) - 1, 16):
                #     int_x = 0
                #     for b in range(7, 0, -1):
                #         hex_to_int = (int(str(hexdata[x + b * 2]), 16)) * 16 + int(str(hexdata[x + b * 2 + 1]), 16)
                #         int_x = (int_x + hex_to_int) << 8
                #
                #     hex_to_int = (int(str(hexdata[x]), 16)) * 16 + int(str(hexdata[x + 1]), 16)  # acr 2017-11-17 this should fix the problem
                #     int_x = (int_x + hex_to_int)

                    if (((int_x & 0xFF00000000000000) >> 59) == 0x00 and self.MODE == 0):
                        mystruct.runNo = self.RUN
                        mystruct.gemroc_id = self.GEMROC_ID
                        mystruct.tiger_id = (int_x>>56)&0x7
                        mystruct.channel_id = (int_x>>48)&0x3F
                        mystruct.tac_id = (int_x>>46)&0x3
                        mystruct.tcoarse = (int_x >> 30)&0xFFFF
                        mystruct.ecoarse = (int_x >> 20)&0x3FF
                        mystruct.tfine = (int_x >> 10)&0x3FF
                        mystruct.efine = int_x & 0x3FF
                        mystruct.tcoarse_10b = (int_x >> 30)&0x3FF


                        if (((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))>0:   
                            mystruct.delta_coarse=(((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF))
                        else:
                            mystruct.delta_coarse=(((int_x >> 20)&0x3FF) - ((int_x >> 30)&0x3FF)) + 1024
                        temp_ecoarse = mystruct.ecoarse
                        mystruct.charge_SH = int_x & 0x3FF

                        if(self.GEMROC_ID<4):
                            mystruct.layer_id = 1
                        elif(self.GEMROC_ID>3):
                            mystruct.layer_id = 2
                        if(self.GEMROC_ID>11):
                            mystruct.layer_id = 0
                        
                        tree.Fill()

                    if(self.MODE == 1):
                        if (((int_x & 0xE000000000000000)>>61) == 0x6):
                            packet_header = 1
                            LOCAL_L1_COUNT_31_6 = int_x >> 32 & 0x3FFFFFF
                            LOCAL_L1_COUNT_5_0  = int_x >> 24 & 0x3F
                            LOCAL_L1_COUNT      = (LOCAL_L1_COUNT_31_6 << 6) + LOCAL_L1_COUNT_5_0
                            LOCAL_L1_TIMESTAMP  = int_x & 0xFFFF
                            pre_pretimestamp = pre_timestamp
                            pre_timestamp = l1timestamp
                            pre_l1count = l1count
                            l1count = LOCAL_L1_COUNT
                            l1timestamp = LOCAL_L1_TIMESTAMP
                            
                        if(((int_x & 0xC000000000000000)>>62) == 0x0 and packet_header == 1 and packet_tailer == 0):
                            LOCAL_L1_TS_minus_TIGER_COARSE_TS = LOCAL_L1_TIMESTAMP - ((int_x >> 32) & 0xFFFF)

                            lstigerid.append((int_x>>59)&0x7)
                            lschannel_id.append((int_x>>50)&0x3F)
                            lstac_id.append((int_x>>48)&0x3)
                            lsecoarse.append((int_x >> 20)&0x3FF)
                            lstfine.append((int_x >> 10)&0x3FF)
                            lsefine.append(int_x & 0x3FF)
                            lslasttigerframenum.append((int_x >> 56)& 0x7)

                            temp_ecoarse = (int_x >> 20)&0x3FF
                            lstcoarse_10b.append(((int_x>>32)&0x3FF))
                            temp_tcoarse = ((int_x>>32)&0x3FF)

                            tcoarse = (int_x>>32)&0xFFFF
                            ecoarse = (int_x>>20)&0x3FF
                            if (((int_x >> 20)&0x3FF) - ((int_x >> 32)&0x3FF))>0:   
                                lsdelta_coarse.append(((int_x >> 20)&0x3FF) - ((int_x >> 32)&0x3FF))
                            else:
                                lsdelta_coarse.append(((int_x >> 20)&0x3FF) - ((int_x >> 32)&0x3FF) + 1024)
                            #lsl1ts_min_tcoarse.append(LOCAL_L1_TIMESTAMP-tcoarse)

                            count_mismatch = 0
                            #if(l1timestamp < 1566):
                                #if(tcoarse > 64150):
                                   # tcoarse = tcoarse - 2**16
                            
                            #if((LOCAL_L1_TIMESTAMP-tcoarse) < 1300 or (LOCAL_L1_TIMESTAMP-tcoarse) > 1600):
                            #    if((pre_timestamp-tcoarse) < 1300 or (pre_timestamp-tcoarse) > 1600):
                            #        count_mismatch = 2
                            #        lsl1ts_min_tcoarse.append(pre_pretimestamp-tcoarse)
                            #        l1count_new.append(l1count-2)
                            #    else:
                            #        count_mismatch = 1
                            #        lsl1ts_min_tcoarse.append(pre_timestamp-tcoarse)
                            #        l1count_new.append(l1count-1)
                            #else:
                            #    lsl1ts_min_tcoarse.append(LOCAL_L1_TIMESTAMP-tcoarse)
                            #    l1count_new.append(l1count)

                            lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse
                            l1count_new_to_append = l1count
                            
                            if(not(((LOCAL_L1_TIMESTAMP - tcoarse) > 1299 and (LOCAL_L1_TIMESTAMP - tcoarse) < 1567) or (LOCAL_L1_TIMESTAMP - tcoarse) < -63970)):
                                flag_swap1 = True
                            if(flag_swap1):
                                if((int_x>>59)&0x7 > 3):
                                    lsl1ts_min_tcoarse_to_append = pre_timestamp - tcoarse
                                    l1count_new_to_append = l1count-1
                                    count_mismatch = 1
                                    if(not(((lsl1ts_min_tcoarse_to_append) > 1299 and (lsl1ts_min_tcoarse_to_append) < 1567) or (lsl1ts_min_tcoarse_to_append) < -63970)): 
                                        flag_swap2 = True
                            
                            if(flag_swap2):
                                if((int_x>>59)&0x7 > 3):
                                    lsl1ts_min_tcoarse_to_append = pre_pretimestamp - tcoarse
                                    l1count_new_to_append = l1count-2
                                    count_mismatch = 2
                                    
                            if(lsl1ts_min_tcoarse_to_append < 0):
                                tcoarse = tcoarse - 2**16
                                if((int_x>>59)&0x7 > 3):
                                    if flag_swap1:
                                        lsl1ts_min_tcoarse_to_append = pre_timestamp -tcoarse
                                    if flag_swap2:
                                        lsl1ts_min_tcoarse_to_append = pre_pretimestamp -tcoarse
                                    if not (flag_swap1 or flag_swap2):
                                        lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse
                                else:
                                    lsl1ts_min_tcoarse_to_append = LOCAL_L1_TIMESTAMP - tcoarse

                            lsl1ts_min_tcoarse.append(lsl1ts_min_tcoarse_to_append)
                            l1count_new.append(l1count_new_to_append)


                            lscount_mismatch.append(count_mismatch)
                            lstcoarse.append(tcoarse)



                            lscharge_SH.append(int_x & 0x3FF)

                        if(((int_x & 0xE000000000000000)>>61) == 0x7):
                            packet_tailer = 1
                            gemrocid = (int_x >> 32)&0x1F
                           

                        #if(((int_x & 0xF000000000000000)>>60) == 0x4):
                            #pre_udp_packet = udp_packet
                            #udp_packet = (((int_x >> 32)&0xFFFFF) + ((int_x >> 0) & 0xFFFFFFF))

                        if(packet_header == 1 and packet_tailer == 1):
                            for x in range(len(lstac_id)):
                                mystruct.channel_id = lschannel_id.pop()
                                mystruct.tac_id = lstac_id.pop()
                                mystruct.tcoarse = lstcoarse.pop()
                                mystruct.ecoarse = lsecoarse.pop()
                                mystruct.tfine = lstfine.pop()
                                mystruct.efine = lsefine.pop()
                                mystruct.tcoarse_10b = lstcoarse_10b.pop()
                                mystruct.charge_SH = lscharge_SH.pop()
                                mystruct.tiger_id = lstigerid.pop()
                                mystruct.l1ts_min_tcoarse = lsl1ts_min_tcoarse.pop()
                                mystruct.lasttigerframenum = lslasttigerframenum.pop()
                                mystruct.count_mismatch = lscount_mismatch.pop()
                                mystruct.count_new = l1count_new.pop()
                                mystruct.delta_coarse = lsdelta_coarse.pop()


                                mystruct.local_l1_count = LOCAL_L1_COUNT
                                mystruct.count_ori = l1count
                                mystruct.count = mystruct.count_new
                                mystruct.timestamp = l1timestamp
                                mystruct.gemroc_id = gemrocid
                                mystruct.runNo = self.RUN
                                if(gemrocid<4):
                                    mystruct.layer_id = 1
                                elif(gemrocid>3):
                                    mystruct.layer_id = 2
                                if(gemrocid>11):
                                    mystruct.layer_id = 0
                                hitcounter = hitcounter + 1
                                if(hitcounter>max_hitcount): 
                                    continue
                                tree.Fill()
                            packet_header = 0
                            packet_tailer = 0


        rootFile.Write()
        rootFile.Close()


import sys
if(len(sys.argv) == 4):
    filename=sys.argv[1]
    gemroc_id=sys.argv[2]
    mode=sys.argv[3]
    run=1
elif(len(sys.argv) == 5):
    filename=sys.argv[1]
    gemroc_id=sys.argv[2]
    mode=sys.argv[3]
    run=sys.argv[4]
else:
    filename=input("Insert file name:")
    gemroc_id=input("Insert Gemroc id:")
    mode=input("MODE:(trigger_less:0 trigger_matched:1)")
    run=1

GEM5=reader(gemroc_id,mode,run)
GEM5.write_root(filename)

print ("Decoding: " + filename)
