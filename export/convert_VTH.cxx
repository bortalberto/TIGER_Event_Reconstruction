#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <string>

void convert_VTH(std::string iname="delta_vth.root", std::string oname="delta_vth_mapped.root"){
    auto mapfile = new TFile("mapping_IHEP_L2_2planari_penta.root");
    auto maptree = (TTree*)mapfile->Get("tree");

    int channel_id, gemroc_id, SW_FEB_id, pos_x, pos_v;

    maptree->SetBranchAddress("channel_id", &channel_id);
    maptree->SetBranchAddress("gemroc_id", &gemroc_id);
    maptree->SetBranchAddress("SW_FEB_id", &SW_FEB_id);
    maptree->SetBranchAddress("pos_x", &pos_x);
    maptree->SetBranchAddress("pos_v", &pos_v);

    int mx[20][8][64];
    int mv[20][8][64];
    memset(mx, -1, sizeof(mx));
    memset(mv, -1, sizeof(mv));

    for (int i = 0; i < maptree->GetEntries(); i++) {
	maptree->GetEntry(i);
	mx[gemroc_id][SW_FEB_id][channel_id] = pos_x;
	mv[gemroc_id][SW_FEB_id][channel_id] = pos_v;
    }
    
    //mock strips for planar GEM
    for (int ch = 0; ch<64; ch++) {
	mx[11][0][ch] = ch;
    mv[11][1][ch] = ch;

    }

    auto datafile = new TFile(iname.c_str());
    auto datatree = (TTree*)datafile->Get("tree");

    int dgemroc, dFEB, dchannel;

    datatree->SetBranchAddress("channel_id", &dchannel);
    datatree->SetBranchAddress("gemroc_id", &dgemroc);
    datatree->SetBranchAddress("software_feb_id", &dFEB);

    auto ofile = new TFile(oname.c_str(),"RECREATE");
    TTree *otree = datatree->CloneTree();

    int strip_x, strip_v;

    auto bstrip_x = otree->Branch("strip_x",&strip_x,"strip_x/I");
    auto bstrip_v = otree->Branch("strip_v",&strip_v,"strip_v/I");

    for (int i = 0; i < datatree->GetEntries(); i++) {
        datatree->GetEntry(i);

        strip_x = mx[dgemroc][dFEB][dchannel];
        strip_v = mv[dgemroc][dFEB][dchannel];
        bstrip_x->Fill();
        bstrip_v->Fill();
    }

    ofile->Write();
    ofile->Close();
    exit(1);
}
