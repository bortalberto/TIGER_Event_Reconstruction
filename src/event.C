#include "event.h"

int count_unique(std::vector<int> v){
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
    return v.size();
}

int count_diff(std::vector<int> v){
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());

    return fabs(v.back() - v.front());
}

void event(){
  std::string iname="ana.root";
  std::string oname="event.root";

    auto file = new TFile(iname.c_str());
    auto tree = (TTree*)file->Get("tree");

    int dchannel, dgemroc, dFEB, dcount, dtimestamp, dstrip_x, dstrip_v, dl1ts_min_tcoarse, dlasttigerframenum, dchip, dFEB_label, drunNo, dlayer, dtrigg_flag;
    float dcharge_SH, dpos_phi, dtcoarse, decoarse, dtfine, define, dttrigg, dtrigg_tcoarse, dconstant, dslope, dqmax, dtime, dradius, ddelta_coarse; 

    tree->SetBranchAddress("runNo",&drunNo);
    tree->SetBranchAddress("layer",&dlayer);
    tree->SetBranchAddress("channel",&dchannel);
    tree->SetBranchAddress("gemroc",&dgemroc);
    tree->SetBranchAddress("FEB",&dFEB);
    tree->SetBranchAddress("charge_SH",&dcharge_SH);
    tree->SetBranchAddress("count",&dcount);
    tree->SetBranchAddress("timestamp",&dtimestamp);
    tree->SetBranchAddress("pos_phi",&dpos_phi);
    tree->SetBranchAddress("radius",&dradius);
    tree->SetBranchAddress("strip_v",&dstrip_v);
    tree->SetBranchAddress("strip_x",&dstrip_x);
    tree->SetBranchAddress("tcoarse",&dtcoarse);
    tree->SetBranchAddress("ecoarse",&decoarse);
    tree->SetBranchAddress("tfine",&dtfine);
    tree->SetBranchAddress("efine",&define);
    tree->SetBranchAddress("l1ts_min_tcoarse",&dl1ts_min_tcoarse);
    tree->SetBranchAddress("trigg_tcoarse",&dtrigg_tcoarse);
    tree->SetBranchAddress("lasttigerframenum",&dlasttigerframenum);
    tree->SetBranchAddress("chip",&dchip);
    tree->SetBranchAddress("FEB_label",&dFEB_label);
    tree->SetBranchAddress("QDCcali_constant",&dconstant);
    tree->SetBranchAddress("QDCcali_slope",&dslope);
    tree->SetBranchAddress("QDCcali_qmax",&dqmax);
    tree->SetBranchAddress("time",&dtime);
    tree->SetBranchAddress("trigg_flag",&dtrigg_flag);
    tree->SetBranchAddress("delta_coarse",&ddelta_coarse);

    std::multimap<int,int> mevt;

    for (int i = 0; i < tree->GetEntries(); i++) {
	tree->GetEntry(i);
	mevt.insert(std::pair<int,int>(dcount,i));
    }

    auto ofile = new TFile(oname.c_str(),"RECREATE");
    auto otree = new TTree("tree","tree");

    int evtNo, nhits, ngemrocs, ntimestamp, runNo, ntcoarse_L1_TP, ntcoarse_L2_TP;
    float trigg_tcoarse;
    int tcount[900], tchannel[900], tgemroc[900], tFEB[900], ttimestamp[900], tstrip_x[900], tstrip_v[900], tl1ts_min_tcoarse[900], tchip[900], tFEB_label[900], tquality[900], tlayer[900];
    float tcharge_SH[900], tpos_phi[900], ttcoarse[900], tecoarse[900], ttfine[900], tefine[900], t_min_ttrigg[900], tconstant[900], tslope[900], tqmax[900], ttime[900], tradius[900], ttrigg[900], delta_coarse[900];

    otree->Branch("runNo",&runNo,"runNo/i"); // the run index
    otree->Branch("evtNo",&evtNo,"evtNo/i"); // the event number index in each run
    otree->Branch("nhits",&nhits,"nhits/i"); // the number of hits in each event
    otree->Branch("ngemrocs",&ngemrocs,"ngemrocs/i"); // the number of gemrocs fired in each event
    otree->Branch("ntimestamp",&ntimestamp,"ntimestamp/i"); // the number of different GEMROC LOCAL_L1_TIMESTAMP values in each event
    otree->Branch("tcoarse_L1_TP_diff",&ntcoarse_L1_TP,"ntcoarse_L1_TP/i"); // the number of different TCOARSEs for test pulse in L1
    otree->Branch("tcoarse_L2_TP_diff",&ntcoarse_L2_TP,"ntcoarse_L2_TP/i"); // the number of different TCOARSEs for test pulse in L2
    otree->Branch("trigg_tcoarse",&trigg_tcoarse,"trigg_tcoarse/f"); // tcoarse value for the trigger channel in each event
    otree->Branch("local_l1_count",&tcount,"tcount[nhits]/i"); // GEMROC LOCAL_L1_COUNT number for each hit
    otree->Branch("layer",&tlayer,"tlayer[nhits]/i"); // Layer No. for each hit
    otree->Branch("channel",&tchannel,"tchannel[nhits]/i");// channel ID for each hit
    otree->Branch("gemroc",&tgemroc,"tgemroc[nhits]/i");// GEMROC ID for each hit
    otree->Branch("FEB_SW",&tFEB,"tFEB[nhits]/i");// FEB ID for each hit
    otree->Branch("local_l1_timestamp",&ttimestamp,"ttimestamp[nhits]/i");// GEMROC LOCAL_L1_TIMESTAMP for each hit
    otree->Branch("charge_SH",&tcharge_SH,"tcharge_SH[nhits]/f");// charge in S&H mode, with QDC calibration, for each hit
    otree->Branch("radius",&tradius,"tradius[nhits]/f");// radius for each hit, L1: 90.223 mm; L2: 129.8 mm
    otree->Branch("pos_phi",&tpos_phi,"tpos_phi[nhits]/f");// phi positon (strip_x) for each hit
    otree->Branch("strip_v",&tstrip_v,"tstrip_v[nhits]/i");// strip v number of each hit
    otree->Branch("strip_x",&tstrip_x,"tstrip_x[nhits]/i");// strip x number of each hit
    otree->Branch("chip",&tchip,"tchip[nhits]/i"); // FEB chip ID for each hit
    otree->Branch("FEB_label",&tFEB_label,"tFEB_label[nhits]/i");//FEB label number for each hit
    otree->Branch("tcoarse",&ttcoarse,"ttcoarse[nhits]/f"); // tcoarse value for each hit
    otree->Branch("ecoarse",&tecoarse,"tecoarse[nhits]/f");// ecoarse value for each hit
    otree->Branch("tfine",&ttfine,"ttfine[nhits]/f");// tfine value for each hit
    otree->Branch("efine",&tefine,"tefine[nhits]/f");// efine value for each hit
    otree->Branch("l1ts_min_tcoarse",&tl1ts_min_tcoarse,"tl1ts_min_tcoarse[nhits]/i"); // GEMROC LOCAL_L1_TIMESTAMP - tcoarse for each hit
    otree->Branch("ttrigg",&ttrigg,"ttrigg[nhits]/f"); // tcoarse of each trigger channel 
    otree->Branch("t_min_ttrigg",&t_min_ttrigg,"t_min_ttrigg[nhits]/f"); // tcoarse of each hit -  tcoarse from trigger channel in each event
    otree->Branch("quality",&tquality,"tquality[nhits]/i"); // tags of good event for each hit
    otree->Branch("QDCcali_constant",&tconstant,"tconstant[nhits]/f"); // constant value from DQC calibration curve for each hit
    otree->Branch("QDCcali_slope",&tslope,"tslope[nhits]/f");// slope value from DQC calibration curve for each hit
    otree->Branch("QDCcali_qmax",&tqmax,"tqmax[nhits]/f");// Q max value from DQC calibration curve for each hit
    otree->Branch("time",&ttime,"ttime[nhits]/f");// Q max value from DQC calibration curve for each hit
    otree->Branch("delta_coarse",&delta_coarse,"delta_coarse[nhits]/f");// Q max value from DQC calibration curve for each hit
    
    std::multimap<int, int>::iterator it = mevt.begin();
    int itcount = it->first;
    nhits = 0;
    evtNo = 0;
    ngemrocs = 0;
    runNo = 0;
    std::vector<int> vgemrocs;
    std::vector<int> vtimestamp;
    std::vector<float> vttrigg;
    std::vector<int> vlasttigerframenum;
    std::vector<int> vtrigg_gemroc;

    std::vector<int> vtcoarse_L1_TP;
    std::vector<int> vtcoarse_L2_TP;

    for (std::pair<int, int> elem : mevt){
	if(itcount==elem.first||nhits==0){ //an event collects all hits which have a same GEMROC LOCAL_COUNT number
	    itcount = elem.first;
	    tree->GetEntry(elem.second);
	    if(!runNo) runNo = drunNo;

	    ttrigg[nhits] = -999;
	    if(dtrigg_flag){ // trigger channgel
		trigg_tcoarse = dtrigg_tcoarse;
		vttrigg.push_back(dtrigg_tcoarse);
		vtrigg_gemroc.push_back(dgemroc);
		vlasttigerframenum.push_back(dlasttigerframenum);
		ttrigg[nhits] = dtrigg_tcoarse;

	    }
	    //if(dcharge_SH<0) continue;
	    //if(dlayer==2&&dchannel>60) continue; // channel 61,62,63 for L2 are not used!
	    //if(dstrip_v==-1&&dstrip_x==-1) continue; // some stripID are not properly setted at the present !!!

	    tcount[nhits] = dcount;
	    tlayer[nhits] = dlayer;
	    tchannel[nhits] = dchannel;
	    tgemroc[nhits] = dgemroc;
	    tFEB[nhits] = dFEB;
	    ttimestamp[nhits] = dtimestamp;
	    tcharge_SH[nhits] = dcharge_SH;
	    tpos_phi[nhits] = dpos_phi;
	    tradius[nhits] = dradius;
	    tstrip_x[nhits] = dstrip_x;
	    tstrip_v[nhits] = dstrip_v;
	    ttcoarse[nhits] = dtcoarse;
	    tecoarse[nhits] = decoarse;
	    ttfine[nhits] = dtfine;
	    tchip[nhits] = dchip;
	    tFEB_label[nhits] = dFEB_label;
	    tefine[nhits] = define;
	    tl1ts_min_tcoarse[nhits] = dl1ts_min_tcoarse;
	    tconstant[nhits] = dconstant;
	    tslope[nhits] = dslope;
	    tqmax[nhits] = dqmax;
	    tquality[nhits] = 0;
	    ttime[nhits] = dtime;
	    delta_coarse[nhits] = ddelta_coarse;
	    ////tag tquality
	    //if(dFEB_label==10&&dchip==2) tquality[nhits] = 1;
	    //if(dFEB_label==13&&dchip==1) tquality[nhits] = 1;
	    //if(dFEB_label==14) tquality[nhits] = 1;
	    //if(dFEB_label==3&&dchip==1) tquality[nhits] = 1;
	    //if(dFEB_label==0&&dchip==1) tquality[nhits] = 1;
	    ////the above tagging only avabile for Ferrara data on Nov. 2018
	    //
	    //check test pulse
	    if(dlayer==1){
		//if(dchannel<36&&dchannel%7==0) vtcoarse_L1_TP.push_back(dtcoarse);
		if(dchannel==20) vtcoarse_L1_TP.push_back(dtcoarse);
	    }
	    if(dlayer==2){
		//if(dchannel==20) vtcoarse_L2_TP.push_back(dtcoarse);
		if(dchannel==20&&dFEB<4&&dgemroc<10) vtcoarse_L2_TP.push_back(dtcoarse);
	    }

	    vgemrocs.push_back(dgemroc);
	    vtimestamp.push_back(dtimestamp);
	    nhits++;
	    if(nhits>900) {std::cout<<"the number of hits is out of ARRAY range"<<std::endl; break;}
	}
	else{
	    ngemrocs = count_unique(vgemrocs);
	    ntimestamp = count_unique(vtimestamp);

	    if(!vtcoarse_L1_TP.empty()) ntcoarse_L1_TP = count_diff(vtcoarse_L1_TP);
	    else ntcoarse_L1_TP = 9999;

	    if(!vtcoarse_L2_TP.empty()) ntcoarse_L2_TP = count_diff(vtcoarse_L2_TP);
	    else ntcoarse_L2_TP = 9999;

	    if(ngemrocs<2 || nhits<2){ //at least two gemrocs should be fired with two hits
		nhits = 0;
		vgemrocs.clear();
		vtimestamp.clear();
		vttrigg.clear();
		vtcoarse_L1_TP.clear();
		vtcoarse_L2_TP.clear();
		continue;
	    }
	    else{
		int good = 0;
		//for(int i=1; i<nhits; i++){
		//  if(sin(tpos_phi[0])*sin(tpos_phi[i])<0) good++; //require one hit with negtive pos_phi, anther one with postive pos_phi
		//}
		for(int i=0; i<nhits; i++){
                    good++; //no pos_phi requirement
                }
		if(good==0){
		    nhits = 0;
		    trigg_tcoarse = -9999;
		    vgemrocs.clear();
		    vtimestamp.clear();
		    vttrigg.clear();
		    continue;
		}

		//for(int i=0; i<vttrigg.size(); i++){
		//    cout<<vtrigg_gemroc[i]<<"\t";
		//}
		//cout<<""<<endl;
		//for(int i=0; i<vttrigg.size(); i++){
		//    cout<<vttrigg[i]<<"\t";
		//}
		//cout<<""<<endl;
		//cout<<""<<endl;


		for(int i=0; i<nhits; i++){
		    //if(tstrip_x[i]<0){ 
		    //    t_min_ttrigg[i]=-9999; 
		    //    continue;
		    //}
		    //if(vttrigg.size()<1) vttrigg.push_back(-99);
		    if(tlayer[i]==1){
			if(ntcoarse_L1_TP>3) t_min_ttrigg[i] = -9999;
			else t_min_ttrigg[i] = TMath::Abs(ttcoarse[i]) - TMath::Abs(vtcoarse_L1_TP[0]);
		    }
		    if(tlayer[i]==2){
			if(ntcoarse_L2_TP>3) t_min_ttrigg[i] = -9999;
			else t_min_ttrigg[i] = TMath::Abs(ttcoarse[i]) - TMath::Abs(vtcoarse_L2_TP[0]);} //if there are two trigger time in one event, only the first one is used here.
		}
	    }

	    evtNo++;
	    otree->Fill();
	    nhits=0;
	    trigg_tcoarse = -9999;
	    vgemrocs.clear();
	    vtimestamp.clear();
	    vttrigg.clear();
	    vtcoarse_L1_TP.clear();
	    vtcoarse_L2_TP.clear();
	    if(evtNo%100==0) std::cout<<"Evt. No.\t"<<evtNo<<std::endl;
	    if(evtNo>1000) break;
	}
	
	}

    otree->Write();
    file->Close();
    ofile->Close();
}
