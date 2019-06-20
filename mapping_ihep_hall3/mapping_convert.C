#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TDirectory.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>

float x_to_phi(int strip_x, int layer){
    float pitch = 0;
    if(layer==1){
	pitch = 2*TMath::Pi()/856.;
    }
    else if(layer==2){
	pitch = 2*TMath::Pi()/1260.;
    }
    return (strip_x - 1)*pitch;
}


//float v_to_phi(int strip_v, int layer){
//    float pitch = 0;
//    if(layer==1){
//        pitch = 2*TMath::Pi()/856.;
//    }
//    else if(layer==2){
//        pitch = 2*TMath::Pi()/1260.;
//    }
//    return (strip_v - 1)*pitch;
//}


void mapping_convert(){
    auto file = new TFile("mapping_IHEP.root","RECREATE");
    auto tree = new TTree("tree","tree");

    int HW_FEB_id,chip_id,ch_id,gemroc_id,SW_FEB_id,layer_id,FEB_label,pos_x,pos_v;
    float phi;
    std::string view;

    tree->Branch("layer_id",&layer_id,"layer_id/I");
    tree->Branch("channel_id",&ch_id,"channel_id/I");
    tree->Branch("gemroc_id",&gemroc_id,"gemroc_id/I");
    tree->Branch("SW_FEB_id",&SW_FEB_id,"SW_FEB_id/I");
    tree->Branch("FEB_label",&FEB_label,"FEB_label/I");
    tree->Branch("HW_FEB_id",&HW_FEB_id,"HW_FEB_id/I");
    tree->Branch("chip_id",&chip_id,"chip_id/I");
    tree->Branch("pos_x",&pos_x,"pos_x/I");
    tree->Branch("pos_v",&pos_v,"pos_v/I");
    tree->Branch("phi",&phi,"phi/F");


    int mpos_x[2][64][44];
    int mpos_v[2][64][44];
    memset(mpos_x, -1, sizeof(mpos_x));
    memset(mpos_v, -1, sizeof(mpos_v));

    for(int layer=0; layer<2; layer++){
	std::string FILENAME = Form("L%d_mapping.txt",layer+1);
	std::ifstream data_file(FILENAME);

	if(!data_file){
	    std::cout << std::string("could not open the file ") + FILENAME<< std::endl;
	    continue;
	}
	//cout<<"!!!!!!"<<endl<<endl<<endl<<endl;
	while(data_file >> chip_id >> ch_id >> view >> FEB_label){
	  if(ch_id==-1) continue;
	  //cout<<chip_id<<" "<<ch_id<<" "<<view<<" "<<FEB_label<<endl;
	  string vista = view.substr(0,1);
	  string svista = view.erase(0,1);
	  //cout<<vista<<" "<<svista<<endl;
	  int channel = std::stoi(svista);
	  if(vista=="x")  mpos_x[chip_id-1][ch_id][FEB_label] = channel;
	    else if(vista=="v") mpos_v[chip_id-1][ch_id][FEB_label] = channel;
	  //cout<<mpos_v[chip_id-1][ch_id][FEB_label]<<endl;

	}
    }


	for(int chip=0; chip<2; chip++){
	    for(ch_id=0; ch_id<64; ch_id++){
		for(FEB_label=0; FEB_label<44; FEB_label++){

		    pos_x = mpos_x[chip][ch_id][FEB_label];
		    pos_v = mpos_v[chip][ch_id][FEB_label];
		    chip_id =  chip + 1;

		    gemroc_id = FEB_label/4;
                    SW_FEB_id = (FEB_label%4)*2 + chip_id-1;

		    if(gemroc_id<4) layer_id = 1;
		    if(gemroc_id>3) layer_id = 2;
		    phi = x_to_phi(pos_x, layer_id);


		    if     (FEB_label==0)  HW_FEB_id = 5;
		    else if(FEB_label==1)  HW_FEB_id = 20;
		    else if(FEB_label==2)  HW_FEB_id = 3; // old: 6;  new:12
		    else if(FEB_label==3)  HW_FEB_id = 13;
		    else if(FEB_label==4)  HW_FEB_id = 14;
		    else if(FEB_label==5)  HW_FEB_id = 3;
		    else if(FEB_label==6)  HW_FEB_id = 18;
		    else if(FEB_label==7)  HW_FEB_id = 19;
		    else if(FEB_label==8)  HW_FEB_id = 24;
		    else if(FEB_label==9)  HW_FEB_id = 11;
		    else if(FEB_label==10) HW_FEB_id = 22;
		    else if(FEB_label==11) HW_FEB_id = 21;
		    else if(FEB_label==12) HW_FEB_id = 17;
		    else if(FEB_label==13) HW_FEB_id = 25;
		    else if(FEB_label==14) HW_FEB_id = 10;
		    else if(FEB_label==15) HW_FEB_id = 23;
		    else if(FEB_label==16) HW_FEB_id = 4;
		    else if(FEB_label==17) HW_FEB_id = 17;
		    else if(FEB_label==18) HW_FEB_id = 18;
		    else if(FEB_label==19) HW_FEB_id = 6;
		    else if(FEB_label==20) HW_FEB_id = 9;
		    else if(FEB_label==21) HW_FEB_id = 33;
		    else if(FEB_label==22) HW_FEB_id = 38;
		    else if(FEB_label==23) HW_FEB_id = 34;
		    else if(FEB_label==24) HW_FEB_id = 37;
		    else if(FEB_label==25) HW_FEB_id = 36;
		    else if(FEB_label==26) HW_FEB_id = 5;
		    else if(FEB_label==27) HW_FEB_id = 1;
		    else if(FEB_label==28) HW_FEB_id = 2;
		    else if(FEB_label==29) HW_FEB_id = 19;
		    else if(FEB_label==30) HW_FEB_id = 52;
		    else if(FEB_label==31) HW_FEB_id = 32; 
		    else if(FEB_label==32) HW_FEB_id = 30;
		    else if(FEB_label==33) HW_FEB_id = 26;
		    else if(FEB_label==34) HW_FEB_id = 27;
		    else if(FEB_label==35) HW_FEB_id = 29;
		    else if(FEB_label==36) HW_FEB_id = 25;
		    else if(FEB_label==37) HW_FEB_id = 11;
		    else if(FEB_label==38) HW_FEB_id = 12;
		    else if(FEB_label==39) HW_FEB_id = 7;
		    else if(FEB_label==40) HW_FEB_id = 20;
		    else if(FEB_label==41) HW_FEB_id = 21;
		    else if(FEB_label==42) HW_FEB_id = 24;
		    else if(FEB_label==43) HW_FEB_id = 3;
		    else                   HW_FEB_id = -1;


		    tree->Fill();
		}
	    }
	}
    file->Write();
    file->Close();

}


