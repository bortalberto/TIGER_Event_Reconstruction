#include <iostream>
#include <fstream>
#include "daq.h"
using namespace std;

int main(int argc, const char* argv[]){
  //string folder = argv[1];
  int run    = atoi(argv[1]);
  int subrun;
  if(argc==3) subrun = atoi(argv[2]);
  if(argc==2) daq(run);
  if(argc==3) daq(run,subrun);
  return 0;
};
