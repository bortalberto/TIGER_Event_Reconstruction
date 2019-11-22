#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TDirectory.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>

void totcalib(){
  auto file = new TFile("ToT_calib.root","RECREATE");
  auto tree = new TTree("tree","tree");

  int vthr;
  float par_a, par_b, par_c, par_d, par_e;

  tree->Branch("vthr",&vthr,"vthr/I");
  tree->Branch("par_a",&par_a,"par_a/F");
  tree->Branch("par_b",&par_b,"par_b/F");
  tree->Branch("par_c",&par_c,"par_c/F");
  tree->Branch("par_d",&par_d,"par_d/F");
  tree->Branch("par_e",&par_e,"par_e/F");

  std::string FILENAME = "ToT_calib.txt";
  std::ifstream data_file(FILENAME);
  if(!data_file){
    std::cout << std::string("could not open the file ") + FILENAME<< std::endl;
    return;
  }
  while(data_file >> vthr >> par_a >> par_b >> par_c >> par_d >> par_e){
    //cout<<vthr<<" "<<par_a<<" "<<par_e<<endl;
    tree->Fill();
  }
  tree->Write();
  file->Close();


  return;
}
