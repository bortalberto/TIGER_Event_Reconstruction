#include <iostream>
#include <fstream>
#include "ext.h"
using namespace std;

int main(int argc, const char* argv[]){
  //string folder = argv[1];
  int run    = atoi(argv[1]);
  int feb, chip, channel;
  if(argc==5) {
    feb     = atoi(argv[2]);
    chip    = atoi(argv[3]);
    channel = atoi(argv[4]);
  }
  if(argc==2) ext(run);
  else if(argc==5) ext(run,feb,chip,channel);
  else cout<<"BAD ARGUMENT"<<endl;
  return 0;
};
