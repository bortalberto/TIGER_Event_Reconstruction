#ifndef ext_h
#define ext_h
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <string>
#include "TCanvas.h"
#include "TH2F.h"
#include "TCut.h"
#include "TChain.h"
#include "TF1.h"
#include <map>
#include "common.h"
#include "TSystem.h"
#include <thread>
#include <chrono>

const string DataDir="/home/ihep_data/data/raw_root/";

using namespace std;
void ext(int);
void ext(int,int,int,int);
void ext_i(int,int,int);
bool init(int);
TString Get_Cut(int,int);
TString Get_Cut(int,int,int);
TString Get_Cut(int,int,int,int);
TString Get_Command(int,TString);
TString Get_Command(int,TString,int);
void extract_noise_and_thr(int,int);
int count_line(int);
#endif
