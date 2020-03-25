#ifndef event_h
#define event_h
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <string>
#include "TCanvas.h"
#include "TH2F.h"
#include "TCut.h"
#include <map>
#include "common.h"

using namespace std;
void event(int,int);
void TP_fill(int,int,int,double);
void TP_test(std::vector<int>,std::vector<int>);
void TP_cout();
void TIGER_count(int,int);
void exit_loop();
int count_unique(std::vector<int> v);
int count_diff(std::vector<int> v);

#endif
