#include <iostream>
#include <fstream>
#include "clean.h"
using namespace std;

int main(int argc, const char* argv[]){
  //string folder = argv[1];
  int run    = atoi(argv[1]);
  clean(run);
  return 0;
};
