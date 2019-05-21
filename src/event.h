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

using namespace std;
void event();
int count_unique(std::vector<int> v);
int count_diff(std::vector<int> v);
#endif
