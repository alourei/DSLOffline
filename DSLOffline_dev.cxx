//Offline routine for TRIUMF DSAM S1582 experiment
//D. Perez Loureiro 2016/05/26
//History: 
//2016/05/27 16:00 Time stamps added
//2016/05/27 18:10 New version with more cuts and spectra
//2016/05/27 21:10 The output file saved with the same name of the run
//2016/05/27 22:00 Added total energy spectra on si form alphas and He3
//2016/05/29 22:00 Added total energy spectra on si form alphas and He3
//2016/05/30 15:00 Added addback spectra
//2016/05/30 19:00 Added addback PID Gamma gated matrix
//2016/05/30 18:09 Progress also shown as a percentage
//2016/05/30 19:09 Added the total Run time;
//2016/05/30 19:09 Added Output Tree;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TApplication.h>
#include <TRint.h>
#include <TCutG.h>
#include <TStopwatch.h>

// Header file for the classes stored in the TTree if any.
#include "/home/grifstor/DSAM/GRSISort/include/TGriffin.h"
#include "/home/grifstor/DSAM/GRSISort/include/TGRSIDetectorHit.h"
#include "/home/grifstor/DSAM/GRSISort/include/TGriffinHit.h"
#include "/home/grifstor/DSAM/GRSISort/include/TSceptar.h"
#include "/home/grifstor/DSAM/GRSISort/include/TSceptarHit.h"
#include "/home/grifstor/DSAM/GRSISort/include/TPaces.h"
#include "/home/grifstor/DSAM/GRSISort/include/TPacesHit.h"
#include "/home/grifstor/DSAM/GRSISort/include/TChannel.h"
#include "/home/grifstor/DSAM/GRSISort/include/TGRSIRunInfo.h"

Int_t main(Int_t argc,Char_t **argv){


  //TApplication theApp("Analysis",0,0);
  TRint theApp("Analysis",0,0);
  TStopwatch *w=new TStopwatch();
  w->Start();
  if(argc<4){
    std::cout<<"Usage: "<< argv[0]<<" <calfile> <cutfile> <input file(s)>"<<std::endl;
    return 0;
  }

  else{
    Bool_t Draw_histos=kFALSE;
    TFile *cut_file=new TFile(argv[2]);
    TCutG *cg=(TCutG*)cut_file->Get("alphas");
    //TCutG *cg2=(TCutG*)cut_file->Get("protons_hiE");
    TCutG *cg2=(TCutG*)cut_file->Get("protons");
    TCutG *cg3=(TCutG*)cut_file->Get("all");
    TCutG *cg4=(TCutG*)cut_file->Get("He3_tot");
    TCutG *cg5=(TCutG*)cut_file->Get("He3_ex");
    if(cut_file&&cg&&cg2&&cg3&&cg4){
      std::cout<<"Found cut "<<cg->GetName()<<" in file "<<cut_file->GetName()<<std::endl;
      std::cout<<"Found cut "<<cg2->GetName()<<" in file "<<cut_file->GetName()<<std::endl;
      std::cout<<"Found cut "<<cg3->GetName()<<" in file "<<cut_file->GetName()<<std::endl;
      std::cout<<"Found cut "<<cg4->GetName()<<" in file "<<cut_file->GetName()<<std::endl;
      std::cout<<"Found cut "<<cg5->GetName()<<" in file "<<cut_file->GetName()<<std::endl;
      cg2->SetLineColor(2);
      cg2->SetLineWidth(2);
      cg->SetLineWidth(2);
      cg2->SetLineColor(2);
      cg3->SetLineColor(6);
      cg3->SetLineWidth(2);
      cg4->SetLineColor(3);
      cg4->SetLineWidth(2);
      cg4->SetLineColor(9);
      cg4->SetLineWidth(2);
    }
    else{
      std::cerr<<"ERROR Cut File or cuts not found!"<<std::endl;
      return -1;
    }
   TChannel *theChannel=new TChannel();
    Int_t nchannels=theChannel->ReadCalFile(argv[1]);
    if( nchannels<=0){
      std::cerr<<"ERROR Calibration File not found!"<<std::endl;
      return -1;
    }
    TChain *ch=new TChain("AnalysisTree");
    Int_t n_files=3;
    Double_t Total_time=0;
    while(n_files<argc){
      TFile f(argv[n_files]);
      TGRSIRunInfo *theRunInfo=(TGRSIRunInfo*)f.Get("TGRSIRunInfo");
      Double_t runtime=theRunInfo->RunLength();
      //Int_t runnumber=theRunInfo->RunNumber();
      Total_time+= runtime;
      f.Close();
      ch->AddFile(argv[n_files]);
      std::cout<<"File "<<argv[n_files]<<" attached"<<std::endl;
      //std::cout<<"File "<<argv[n_files]<<" is "<<runnumber<<" and is "<<runtime<<" s long"<<std::endl;
      n_files++;
    }
    //Create the Objects to hold the branches
    TGriffin *grif=new TGriffin();
    TSceptar *scep=new TSceptar();
    TPaces *pace=new TPaces();

    ch->SetBranchAddress("TGriffin",&grif);
    ch->SetBranchAddress("TSceptar",&scep);
    ch->SetBranchAddress("TPaces",&pace);


    Char_t outname[256];
    sprintf(outname,"Test_%s",argv[n_files-1]);
    std::cout<<"Histograms will be saved in "<<outname<<std::endl;
    
    TFile *output_file = new TFile(outname,"RECREATE");
 
    //Histograms of interest
    output_file->cd();
    TDirectory *Histo_dir=output_file->mkdir("Histograms");
    Histo_dir->cd();
    //PID
    TH2F *PID=new TH2F("PID","#DeltaE-E PID plot;E;#Delta E",2500,0,5000,2500,0,5000);
    TH2F *PID_gated=new TH2F("PID_g","#DeltaE-E PID plot gated on 2810 keV;E;#Delta E",200,0,5000,200,0,5000);
    TH2F *PID_bkg=new TH2F("PID_b","#DeltaE-E PID plot  gated on background;E;#Delta E",200,0,5000,200,0,5000);
    //Griffin Clover Detectors
    TH1F *hSiE=new TH1F("hSiE","Energy  loss in 1 mm Silicon",8000,0,8000);
    TH1I *hSiE_multiplicity=new TH1I("hSiEM","Multiplicity in 1 mm Silicon",15,0,15);
    TH1I *hGriffin_multiplicity=new TH1I("hGriffM","Multiplicity in Griffin",15,0,15);
    TH1I *hSiDeltaE_multiplicity=new TH1I("hSiDEM","Multiplicity in 500 #um Silicon",15,0,15);
    TH2I *hGrifM_SilM=new TH2I("hGriff_SilM","Multiplicity in Griffin VS silicon;Griffin Mult.;Silicon Mult.",15,0,15,15,0,15);
    TH1F *hSiDeltaE=new TH1F("hSiDeltaE","Energy loss in 0.5 mm Silicon",8000,0,8000);
    TH1F *hSiSum_alpha=new TH1F("hSiSum_a","Total Energy in Silicon Detectors gated on #alpha particles",1000,0,7000);
    TH1F *hSiSum_He3=new TH1F("hSiSum_h","Total Energy in Silicon Detectors gated on ^{3}He",500,0,7000);
    TH1F *hGriffin[8];
    TH1F *hGriffin_time[8];
    TH2F *hGriffin_Energy_time[8];
    TH1F *hGriffing_alpha[8];
    TH1F *hGriffing_proton[8];
    TH1F *hGriffing_He3_all[8];
    TH1F *hGriffing_He3_ex[8];
    TH1F *hGriffing_all[8];
    TH1F *hGriffing_anti[8];
    //Total Spectra
    TH1F *hGriffin0=new TH1F("hGriff_0","Total Energy in Griffin 0 deg",8000,0,8000);
    TH1F *hGriffin90=new TH1F("hGriff_90","Total Energy in Griffin 90 deg",8000,0,8000);
    TH1F *hGriffin0g_alpha=new TH1F("hGriff_0g_alphas","Total Energy in Griffin 0 deg Gated on #alpha particles",8000/4,0,8000);
    TH1F *hGriffin90g_alpha=new TH1F("hGriff_90g_alphas","Total Energy in Griffin 90 deg Gated on #alpha particles",8000/4,0,8000);
    TH1F *hGriffin0g_proton=new TH1F("hGriff_0g_protons","Total Energy in Griffin 0 deg Gated on protons",8000/4,0,8000);
    TH1F *hGriffin90g_proton=new TH1F("hGriff_90g_protons","Total Energy in Griffin 90 deg Gated on protons",8000/4,0,8000);
    TH1F *hGriffin0g_all=new TH1F("hGriff_0g_all","Total Energy in Griffin 0 deg Gated on PID",8000/4,0,8000);
    TH1F *hGriffin90g_all=new TH1F("hGriff_90g_all","Total Energy in Griffin 90 deg Gated on PID",8000/4,0,8000);
    TH1F *hGriffin0g_He3_all=new TH1F("hGriff_0g_He3_all","Total Energy in Griffin 0 deg Gated on ^{3}He",8000/4,0,8000);
    TH1F *hGriffin90g_He3_all=new TH1F("hGriff_90g_He3_all","Total Energy in Griffin 90 deg Gated on ^{3}He",8000/4,0,8000);
    TH1F *hGriffin0g_He3_ex=new TH1F("hGriff_0g_He3_ex","Total Energy in Griffin 0 deg Gated on ^{3}He 1st excited state",8000/4,0,8000);
    TH1F *hGriffin90g_He3_ex=new TH1F("hGriff_90g_He3_ex","Total Energy in Griffin 90 deg Gated on ^{3}He 1st excited state",8000/4,0,8000);
    TH1F *hGriffin0g_anti=new TH1F("hGriff_0g_anti","Total Energy in Griffin 0 deg Anticoincidence with PID",8000,0,8000);
    TH1F *hGriffin90g_anti=new TH1F("hGriff_90g_anti","Total Energy in Griffin 90 deg Anticoincidence with PID",8000,0,8000);
    hGriffin0g_alpha->SetLineColor(1);
    hGriffin90g_alpha->SetLineColor(1);
    hGriffin0g_proton->SetLineColor(2);
    hGriffin90g_proton->SetLineColor(2);
    hGriffin0g_all->SetLineColor(6);
    hGriffin90g_all->SetLineColor(6);
    hGriffin0g_He3_all->SetLineColor(9);
    hGriffin90g_He3_all->SetLineColor(9);
    hGriffin0g_He3_ex->SetLineColor(8);
    hGriffin90g_He3_ex->SetLineColor(8);
    Char_t hname[256];
    Char_t htitle[256];
    for(Int_t j=0;j<8;j++){
      sprintf(hname,"hGriffin_%d",j);
      sprintf(htitle,"Energy in Griffin Germanium detector #%d",j);
      hGriffin[j]=new TH1F(hname,htitle,8000,0,8000);
      sprintf(hname,"hGriffing_alphas%d",j);
      sprintf(htitle,"Energy in Griffin Germanium detector #%d gated on #alpha particles",j);
      hGriffing_alpha[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffing_alpha[j]->SetLineColor(1);
      sprintf(hname,"hGriffing_protons%d",j);
      sprintf(htitle,"Energy in Griffin Germanium detector #%d gated on protons",j);
      hGriffing_proton[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffing_proton[j]->SetLineColor(2);
      sprintf(hname,"hGriffing_all%d",j);
      sprintf(htitle,"Energy in Griffin Germanium detector #%d gated on PID",j);
      hGriffing_all[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffing_all[j]->SetLineColor(6);
      sprintf(hname,"hGriffing_He3_all%d",j);
      sprintf(htitle,"Energy in Griffin Germanium detector #%d gated on ^{3}He",j);
      hGriffing_He3_all[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffing_He3_all[j]->SetLineColor(3);
     sprintf(hname,"hGriffing_He3_ex%d",j);
     sprintf(htitle,"Energy in Griffin Germanium detector #%d gated on ^{3}He 1st excited estate",j);
     hGriffing_He3_ex[j]=new TH1F(hname,htitle,8000/4,0,8000);
     hGriffing_He3_ex[j]->SetLineColor(3);
     sprintf(hname,"hGriffing_times%d",j);
     sprintf(htitle,"Time difference in Griffin Germanium detector with 1mm Si #%d",j);
     hGriffin_time[j]=new TH1F(hname,htitle,1500,-1500,1500);
     hGriffin_time[j]->SetLineColor(4);
     sprintf(hname,"hGriffing_anti%d",j);
     sprintf(htitle,"Energy in Griffin  anticoincidence with PID #%d",j);
     hGriffing_anti[j]=new TH1F(hname,htitle,8000,0,8000);
     hGriffing_anti[j]->SetLineColor(12);
     sprintf(hname,"hGriffin_E_times%d",j);
     sprintf(htitle,"Energy vs Time difference in Griffin #%d",j);
     hGriffin_Energy_time[j]=new TH2F(hname,htitle,1500,-1500,1500,2000,0,8000);
    }

    //Adding Addback spectra
    TH1F *hGriffinAddback[2];
    TH2F *hGriffinAddback_time[2];
    TH1F *hGriffinAddback_alphas[2];
    TH2F *hGriffinAddback_time_alphas[2];
    TH1F *hGriffinAddback_protons[2];
    TH2F *hGriffinAddback_time_protons[2];
    for(Int_t j=0;j<2;j++){
      sprintf(hname,"hGriffinAddback_%d",j);
      sprintf(htitle,"Addback Energy Detector #%d",j+1);
      hGriffinAddback[j]=new TH1F(hname,htitle,8000,0,8000);
      hGriffinAddback[j]->SetLineColor(2);
      hGriffinAddback[j]->SetFillColor(6);
      sprintf(hname,"hGriffinAddback_time%d",j);
      sprintf(htitle,"Addback Energy Detector vs timestamp #%d",j+1);
      hGriffinAddback_time[j]=new TH2F(hname,htitle,1500,-1500,1500,4000,0,8000);
      
      sprintf(hname,"hGriffinAddback_alpha%d",j);
      sprintf(htitle,"Addback Energy Detector #%d gated on #alpha particles",j+1);
      hGriffinAddback_alphas[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffinAddback_alphas[j]->SetLineColor(2);
      hGriffinAddback_alphas[j]->SetFillColor(3);
      sprintf(hname,"hGriffinAddback_time_alphas%d",j);
      sprintf(htitle,"Addback Energy Detector #%d vs timestamp gated on #alpha particles",j+1);
      hGriffinAddback_time_alphas[j]=new TH2F(hname,htitle,1500,-1500,1500,4000/2,0,8000);
      
      sprintf(hname,"hGriffinAddback_proton%d",j);
      sprintf(htitle,"Addback Energy Detector #%d gated on protons",j+1);
      hGriffinAddback_protons[j]=new TH1F(hname,htitle,8000/4,0,8000);
      hGriffinAddback_protons[j]->SetLineColor(2);
      hGriffinAddback_protons[j]->SetFillColor(4);
      sprintf(hname,"hGriffinAddback_time_protons%d",j);
      sprintf(htitle,"Addback Energy Detector #%d vs timestamp gated on protons",j+1);
      hGriffinAddback_time_protons[j]=new TH2F(hname,htitle,1500,-1500,1500,4000/2,0,8000);
      

    }

    //Variables of interest:
    //Griffin
    const Int_t MAX_MULT=20;
    Double_t Energy_Griffin[MAX_MULT];
    UChar_t Channel_Griffin[MAX_MULT];
    ULong_t Timestamp_Griffin[MAX_MULT];
    UChar_t mult_Griffin;
    //Addback
    Double_t Energy_Addback[MAX_MULT];
    UChar_t Detector_Addback[MAX_MULT];
    ULong_t Timestamp_Addback[MAX_MULT];
    UChar_t mult_Addback;
    //DeltaE Silicon
    Double_t Energy_DeltaESil[MAX_MULT];
    ULong_t Timestamp_DeltaESil[MAX_MULT];
    UChar_t mult_DeltaESil;
    //E Silicon
    Double_t Energy_ESil[MAX_MULT];
    ULong_t Timestamp_ESil[MAX_MULT];
    UChar_t mult_ESil;
    
    //Adding a Tree to hold the Events
    output_file->cd();
    //Histo_dir=output_file->mkdir("Tree");
    //Histo_dir->cd();
    TTree *theTree=new TTree("DSLTree","The DSL Events Tree");
    theTree->Branch("Mult_Griffin",&mult_Griffin,"Mult_Griffin/b");
    theTree->Branch("Energy_Griffin",&Energy_Griffin,"EG[Mult_Griffin]/D");
    theTree->Branch("Channel_Griffin",&Channel_Griffin,"CG[Mult_Griffin]/b");
    theTree->Branch("Timestamp_Griffin",&Timestamp_Griffin,"tG[Mult_Griffin]/l");
    
    theTree->Branch("Mult_Addback",&mult_Addback,"Mult_Addback/b");
    theTree->Branch("Energy_Addback",&Energy_Addback,"EA[Mult_Addback]/D");
    theTree->Branch("Detector_Addback",&Detector_Addback,"CA[Mult_Addback]/b");
    theTree->Branch("Timestamp_Addback",&Timestamp_Addback,"tA[Mult_Addback]/l");
    
    theTree->Branch("Mult_DeltaESil",&mult_DeltaESil,"Mult_DeltaESil/b");
    theTree->Branch("Energy_DeltaESil",&Energy_DeltaESil,"ES[Mult_DeltaESil]/D");
    theTree->Branch("Timestamp_DeltaESil",&Timestamp_DeltaESil,"tS[Mult_DeltaESil]/l");
    
    theTree->Branch("Mult_ESil",&mult_ESil,"Mult_ESil/b");
    theTree->Branch("Energy_ESil",&Energy_ESil,"ED[Mult_ESil]/D");
    theTree->Branch("Timestamp_ESil",&Timestamp_ESil,"tD[Mult_ESil]/l");
    
    Long64_t nentries = ch->GetEntries();
    Int_t fiveper=(Int_t)(nentries*0.05);
    //Long64_t nentries =200;
    Long64_t nbytes = 0, nb = 0;
    Int_t n_alphas=0;
    Int_t n_protons=0;
    Int_t n_reactions=0;
    Int_t n_He3=0;
    Int_t n_He3_ex=0;
    //Here is the Event loop
    //Double_t SiliconE=0;  
    //Double_t SilicondE=0;
    //Tracking high multiplicity events in Si
    Int_t multihitThinSi=0;
    Int_t multihitThickSi=0;
    Int_t singlehitThinSi=0;
    Int_t singlehitThickSi=0;
    Int_t Multiplicity=2;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Double_t SiliconE=0;  
      ULong_t Silicon_time=0;  
      Double_t SilicondE=0;  
      Int_t SilMultiplicity=0;
      nb = ch->GetEntry(jentry);   nbytes += nb;
      if(jentry!=0&&jentry%fiveper==0)
	std::cout<<jentry<<"/"<<nentries<<" "<<(Int_t)(jentry/fiveper)*5<<"% of Events processed \r"<<std::flush;
      Int_t nHits=pace->GetMultiplicity();
      mult_DeltaESil=nHits;
      hSiDeltaE_multiplicity->Fill(nHits);
      ULong_t *ThinSilicon_timestamp=new ULong_t[nHits];
      Double_t *ThinSiliconEnergy=new Double_t[nHits];  
      //Double_t ThresholddE=20;  
      if(nHits>0){
	//nHits=1;
	if(nHits==1)
	  singlehitThinSi++;
	if(nHits>=Multiplicity)
	  multihitThinSi++;
	for(Int_t i=0;i<nHits;i++){
	  //std::cout<<nHits<<" in Thin Silicon"<<std::endl;
	  if(pace->GetPacesHit(i)->GetEnergy()>5 && pace->GetPacesHit(i)->GetAddress()==12){
	    ThinSiliconEnergy[i]=pace->GetPacesHit(i)->GetEnergy();
	    Energy_DeltaESil[i]=ThinSiliconEnergy[i];
	    ThinSilicon_timestamp[i]=pace->GetPacesHit(i)->GetTimeStamp();
	    Timestamp_DeltaESil[i]=ThinSilicon_timestamp[i];
	    if(ThinSiliconEnergy[i]>SilicondE)
	      SilicondE=ThinSiliconEnergy[i];
	  }
	}
     
	if(SilicondE>0)
	  hSiDeltaE->Fill(SilicondE);
	
      }
      //std::cout<<"TimeStamp here is "<<The_timeStamp<<std::endl;
      nHits=scep->GetMultiplicity();
      SilMultiplicity=nHits;
      ULong_t *ThickSilicon_timestamp=new ULong_t[nHits];
      Double_t *ThickSiliconEnergy=new Double_t[nHits];  
      //Double_t ThresholdE=80;  
      hSiE_multiplicity->Fill(nHits);
      mult_ESil=nHits;
      
      if(nHits>0){
        if(nHits==1)
	  singlehitThickSi++;
        if(nHits>=Multiplicity)
	  multihitThickSi++;
	//nHits=1;
	for(Int_t i=0;i< nHits;i++){
	  //std::cout<<nHits<<" in Thick Silicon"<<std::endl;
	  if(scep->GetSceptarHit(i)->GetEnergy()>5 && scep->GetSceptarHit(i)->GetAddress()==14){
	    ThickSiliconEnergy[i]=scep->GetSceptarHit(i)->GetEnergy();
	    Energy_ESil[i]=ThickSiliconEnergy[i];
	    ThickSilicon_timestamp[i]=scep->GetSceptarHit(i)->GetTimeStamp();
	    Timestamp_ESil[i]=ThickSilicon_timestamp[i];
	    if(ThickSiliconEnergy[i]>SiliconE){
	      SiliconE=ThickSiliconEnergy[i];
	      Silicon_time=ThickSilicon_timestamp[i];
	    }
	  }
	}
	
	if(SiliconE>0)
	  hSiE->Fill(SiliconE);	 
      }    
      //std::cout<<"TimeStamp here is "<<The_timeStamp<<std::endl;
      //if(ThickSiliconEnergy[i]>5&&ThinSiliconEnergy[i]>5){
      if(SiliconE>5&&SilicondE>5){
	//hSiDeltaE->Fill(ThinSiliconEnergy);
	//hSiE->Fill(ThickSiliconEnergy);
	//PID->Fill(ThickSiliconEnergy[i],ThinSiliconEnergy[i]);
	PID->Fill(SiliconE,SilicondE);
      }
      
      
      nHits=grif->GetMultiplicity();
      mult_Griffin=nHits;
      hGriffin_multiplicity->Fill(nHits);
      hGrifM_SilM->Fill(nHits, SilMultiplicity);
      ULong_t *Griffin_timestamp=new ULong_t[nHits];
      //std::cout<<"Number of hits "<<nHits<<std::endl;
      for(Int_t i=0;i< nHits;i++){
	Double_t  Energy=grif->GetGriffinHit(i)->GetEnergy();
	//Energy_Griffin.push_back(Energy);
	Energy_Griffin[i]=Energy;
	Griffin_timestamp[i]=grif->GetGriffinHit(i)->GetTimeStamp();
	//Timestamp_Griffin.push_back(Griffin_timestamp[i]);
	Timestamp_Griffin[i]=Griffin_timestamp[i];
	Long_t timeStampDifference=Griffin_timestamp[i]-Silicon_time;
	Int_t address=grif->GetGriffinHit(i)->GetAddress();
	//Channel_Griffin.push_back(address);
	Channel_Griffin[i]=address;
	//std::cout<<"Energy "<<Energy<<" Address "<<address<<std::endl;
	hGriffin[address]->Fill(Energy);
	hGriffin_time[address]->Fill(timeStampDifference);
	hGriffin_Energy_time[address]->Fill(timeStampDifference,Energy);
	if(address<4)
	  hGriffin0->Fill(Energy);
	else
	  hGriffin90->Fill(Energy);
	
	//if(timeStampDifference<20){
	if(Energy>0){
	  //if(cg->IsInside(SiliconE,SilicondE)&&SiliconE>3000){
	  if(cg->IsInside(SiliconE,SilicondE)){
	    //	 std::cout<<"Energy in  Silicon "<<ThickSiliconEnergy[0]<<std::endl;
	    //if(ThickSiliconEnergy[0]>1||ThickSiliconEnergy[1]>1){
	    n_alphas++;
	    hSiSum_alpha->Fill(SiliconE+SilicondE);
	    //Double_t sum=(ThickSiliconEnergy[0]+ThinSiliconEnergy[0]);
	    // if(sum>3400&&sum<3600){
	    hGriffing_alpha[address]->Fill(Energy);
	    if(address<4)
	      hGriffin0g_alpha->Fill(Energy);
	    else
	      hGriffin90g_alpha->Fill(Energy);
	  }
	  //}
	  
	  if(cg2->IsInside(SiliconE,SilicondE)){
	    n_protons++;
	    hGriffing_proton[address]->Fill(Energy);
	    if(address<4)
	      hGriffin0g_proton->Fill(Energy);
	    else
	      hGriffin90g_proton->Fill(Energy);
	  }
	  
	  if(cg3->IsInside(SiliconE,SilicondE)){
	    n_reactions++;
	    hGriffing_all[address]->Fill(Energy);
	    if(address<4)
	      hGriffin0g_all->Fill(Energy);
	    else
	      hGriffin90g_all->Fill(Energy);
	  }
	  
	  if(cg4->IsInside(SiliconE,SilicondE)){
	    n_He3++;
	    hSiSum_He3->Fill(SiliconE,SilicondE);
	    hGriffing_He3_all[address]->Fill(Energy);
	    if(address<4)
	      hGriffin0g_He3_all->Fill(Energy);
	    else
	      hGriffin90g_He3_all->Fill(Energy);
	  }
	  
	  if(cg5->IsInside(SiliconE,SilicondE)){
	    n_He3_ex++;
	    hGriffing_He3_ex[address]->Fill(Energy);
	    if(address<4)
	      hGriffin0g_He3_ex->Fill(Energy);
	    else
	      hGriffin90g_He3_ex->Fill(Energy);
	  }
	}//Energy
	//}//TimeStamp
	
       if(!cg3->IsInside(SiliconE,SilicondE)){
	 hGriffing_anti[address]->Fill(Energy);
	 if(address<4)
	   hGriffin0g_anti->Fill(Energy);
	 else
	   hGriffin90g_anti->Fill(Energy);
       }
  
      }
      
      //Adding Addback
      Int_t addbackHits=grif->GetAddbackMultiplicity();
      mult_Addback=addbackHits;
      if(nHits>0){
	//std::cout<<"Number of  Addback hits "<<nHits<<std::endl;
	for(Int_t i=0;i< addbackHits;i++){
	  Double_t Energy=grif->GetAddbackHit(i)->GetEnergy();
	  Energy_Addback[i]=Energy;
	  UInt_t Detector=grif->GetAddbackHit(i)->GetDetector();
	  Detector_Addback[i]=Detector;
	  ULong_t Addback_timestamp=grif->GetAddbackHit(i)->GetTimeStamp();
	  Timestamp_Addback[i]=Addback_timestamp;
	  
	  Long_t difftime=Addback_timestamp-Silicon_time;
	  //std::cout<<Detector<<" "<<Energy<<" "<<Silicon_time<<" "<<difftime<<std::endl;
	  if(Energy>0){
	    hGriffinAddback[Detector-1]->Fill(Energy);
	    hGriffinAddback_time[Detector-1]->Fill(difftime,Energy);
	    
	    //if(cg->IsInside(SiliconE,SilicondE)&&SiliconE>3000){
	    if(cg->IsInside(SiliconE,SilicondE)){
	      hGriffinAddback_alphas[Detector-1]->Fill(Energy);
	      hGriffinAddback_time_alphas[Detector-1]->Fill(difftime,Energy);	      
	    }

	    if(cg2->IsInside(SiliconE,SilicondE)){
	      hGriffinAddback_protons[Detector-1]->Fill(Energy);
	      hGriffinAddback_time_protons[Detector-1]->Fill(difftime,Energy);
	    }
	    
	    if(Detector==1&&Energy>2808&&Energy<2820&&SiliconE>0&&SilicondE>0)
	      PID_gated->Fill(SiliconE,SilicondE);
	    else if(Detector==1&&Energy>2788&&Energy<2808&&SiliconE>0&&SilicondE>0)
	      PID_bkg->Fill(SiliconE,SilicondE);
	    
	  }
	//std::cout<<Detector<<" "<<Energy<<std::endl;
	}
      }

      delete [] ThinSilicon_timestamp;
      delete [] ThinSiliconEnergy;
      delete [] ThickSilicon_timestamp;
      delete [] ThickSiliconEnergy;
      delete [] Griffin_timestamp;
      //std::cout<<"TimeStamp here is "<<The_timeStamp<<std::endl;
      
      theTree->Fill();      
    }//End of event Loop
    
    //output_file->cd();
    //theTree->Print();
    //theTree->Write("",TObject::kOverwrite);
    
    if(Draw_histos){ 
      // TCanvas *c1=new TCanvas("c1","Singles Germanium");
      // c1->Divide(2,4);
      // for(Int_t i=0;i<8;i++){
      //   c1->cd(i+1);
      //   hGriffin[i]->Draw();
      // }
      TCanvas *c2=new TCanvas("c2","PID");
      c2->cd();
      PID->Draw("colz");
      PID->GetXaxis()->SetRangeUser(0,5000);
      PID->GetYaxis()->SetRangeUser(0,4500);
      cg->Draw("same");
      cg2->Draw("same");
      cg3->Draw("same");
      cg4->Draw("same");
      cg5->Draw("same");
      std::cout<<std::endl;
      
      
      
      // TCanvas *c4=new TCanvas("c4","Silicons");
      // c4->cd();
      // c4->Divide(1,2);
      // c4->cd(1);
      // hSiDeltaE->Draw();
      // c4->cd(2);
      // hSiE->Draw();
      
      // TCanvas *c3=new TCanvas("c3","Gated Germanium");
      // c3->Divide(2,4);
      // for(Int_t i=0;i<8;i++){
      //   c3->cd(i+1);
      //   hGriffing_alpha[i]->Draw();
      // }
      
      // TCanvas *c5=new TCanvas("c5","Total Griffin");
      // c5->cd();
      // c5->Divide(1,2);
      // c5->cd(1);
      // hGriffin0->Draw();
      // c5->cd(2);
      // hGriffin90->Draw();
      
      TCanvas *c6=new TCanvas("c6","Total Griffin Gated on alphas");
      c6->cd();
      c6->Divide(1,2);
      c6->cd(1);
      hGriffin0g_alpha->Draw();
      c6->cd(2);
      hGriffin90g_alpha->Draw();
      
      TCanvas *c7=new TCanvas("c7","Total Griffin Gated on protons");
      c7->cd();
      c7->Divide(1,2);
      c7->cd(1);
      hGriffin0g_proton->Draw();
      c7->cd(2);
      hGriffin90g_proton->Draw();
      
      // TCanvas *c8=new TCanvas("c8","Total Griffin Gated on PID");
      // c8->cd();
      // c8->Divide(1,2);
      // c8->cd(1);
      // hGriffin0g_all->Draw();
      // c8->cd(2);
      // hGriffin90g_all->Draw();
      
      // TCanvas *c9=new TCanvas("c9","Total energy deposited in Silicon Detectors #alphas");
      // c9->cd();
      // hSiSum_alpha->Draw();
      
      // TCanvas *c10=new TCanvas("c10","Total Griffin Gated on 3He");
      // c10->cd();
      // c10->Divide(1,2);
      // c10->cd(1);
      // hGriffin0g_He3_all->Draw();
      // c10->cd(2);
      // hGriffin90g_He3_all->Draw();
      
      // TCanvas *c11=new TCanvas("c11","Total Griffin Gated on 3He 1st excited state");
      // c11->cd();
      // c11->Divide(1,2);
      // c11->cd(1);
      // hGriffin0g_He3_ex->Draw();
      // c11->cd(2);
      // hGriffin90g_He3_ex->Draw();
      
      // TCanvas *c12=new TCanvas("c12","Total energy deposited in Silicon Detectors gated on ^{3}He");
      // c12->cd();
      // hSiSum_He3->Draw();
      
      // TCanvas *c13=new TCanvas("c13","TimeStamps Germanium");
      // c13->Divide(2,4);
      // for(Int_t i=0;i<8;i++){
      //   c13->cd(i+1);
      //   hGriffin_time[i]->Draw();
      // }
      
      // TCanvas *c14=new TCanvas("c14","Total Griffin Anticoincidence");
      // c14->cd();
      // c14->Divide(1,2);
      // c14->cd(1);
      // hGriffin0g_anti->Draw();
      // c14->cd(2);
      // hGriffin90g_anti->Draw();
      
      
      // TCanvas *c15=new TCanvas("c15","TimeStamps and Energies Germanium ");
      // c15->Divide(2,4);
      // for(Int_t i=0;i<8;i++){
      //   c15->cd(i+1);
      //   hGriffin_Energy_time[i]->Draw("colz");
      // }
      
      
      TCanvas *c16=new TCanvas("c16","Addback Energies Germanium ");
      c16->Divide(1,2);
      for(Int_t i=0;i<2;i++){
	c16->cd(i+1);
	hGriffinAddback[i]->Draw();
      }
      
      // TCanvas *c17=new TCanvas("c17","Addback Energies vs time Germanium ");
      // c17->Divide(1,2);
      // for(Int_t i=0;i<2;i++){
      //   c17->cd(i+1);
      //   hGriffinAddback_time[i]->Draw("colz");
      // }
      
      
      TCanvas *c18=new TCanvas("c18","Addback Energies Germanium gated alphas ");
      c18->Divide(1,2);
      for(Int_t i=0;i<2;i++){
	c18->cd(i+1);
	hGriffinAddback_alphas[i]->Draw();
      }
      
      // TCanvas *c19=new TCanvas("c19","Addback Energies vs time Germanium gated alphas ");
      // c19->Divide(1,2);
      // for(Int_t i=0;i<2;i++){
      //   c19->cd(i+1);
      //   hGriffinAddback_time_alphas[i]->Draw("colz");
      // }
      
      TCanvas *c20=new TCanvas("c20","Addback Energies Germanium gated protons ");
      c20->Divide(1,2);
      for(Int_t i=0;i<2;i++){
	c20->cd(i+1);
	hGriffinAddback_protons[i]->Draw();
      }
      
      // TCanvas *c21=new TCanvas("c21","Addback Energies vs time Germanium gated protons ");
      // c21->Divide(1,2);
      // for(Int_t i=0;i<2;i++){
      //   c21->cd(i+1);
      //   hGriffinAddback_time_protons[i]->Draw("colz");
      // }
      
      
      TCanvas *c22=new TCanvas("c21"," PID Gated on gammas ");
      c22->Divide(1,2);
      c22->cd(1);
      PID_gated->Draw("colz");
      c22->cd(2);
      PID_bkg->Draw("colz");
      
      
    }//Draw_histos

    std::cout<<"Alphas inside the Gate "<<cg->IntegralHist(PID)<<std::endl;
    std::cout<<"Protons inside the Gate "<<cg2->IntegralHist(PID)<<std::endl;
    std::cout<<"Reactions inside the Gate "<<cg3->IntegralHist(PID)<<std::endl;
    std::cout<<"^{3}He inside the Gate "<<cg4->IntegralHist(PID)<<std::endl;
    std::cout<<"^{3}He_ex inside the Gate "<<cg5->IntegralHist(PID)<<std::endl;
    

    
    std::cout<<"Multiplicity "<<Multiplicity<<" Events in thin Silicon "<<multihitThinSi<<std::endl;
    std::cout<<"Multiplicity "<<Multiplicity<<" Events in thick Silicon "<<multihitThickSi<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Singlehit Events in thin Silicon "<<singlehitThinSi<<std::endl;
    std::cout<<"Singlehit Events in thick Silicon "<<singlehitThickSi<<std::endl;
    
    // new TCanvas();
    // hGriffin_multiplicity->SetLineColor(3);
    // hGriffin_multiplicity->Draw("same");
    // hSiE_multiplicity->Draw();
    // hSiDeltaE_multiplicity->SetLineColor(2);
    // hSiDeltaE_multiplicity->Draw("same");
    // hGrifM_SilM->Draw("colz");
    Int_t hours=Total_time/3600;
    Int_t minutes=(Total_time-hours*3600)/60;
    Int_t seconds=(Total_time-hours*3600-minutes*60);
    
    //std::cout<<"The duration of the run is "<<Total_time<<" seconds"<<std::endl;
    std::cout<<"The duration of the run is "<<hours<<"h "<<minutes<<"min "<<seconds<<"s"<<std::endl;
    
    //theTree->Write();
    output_file->Write("",TObject::kOverwrite);
    w->Stop();
    Double_t realtime=w->RealTime();
    Double_t cputime=w->CpuTime();
    if(realtime<60)
      std::cout<<"The processing time is "<<realtime<<" s"<<std::endl;
    else{
      Int_t minutes=realtime/60;
      Int_t seconds=(realtime-minutes*60);
      std::cout<<"The processing time is "<<minutes<<"min "<<seconds<<" s"<<std::endl;
    }

      TCanvas *co=new TCanvas("c2","PID");
      co->cd();
      PID->Draw("colz");
      PID->GetXaxis()->SetRangeUser(0,5000);
      PID->GetYaxis()->SetRangeUser(0,4500);
      cg->Draw("same");
      cg2->Draw("same");
      cg3->Draw("same");
      cg4->Draw("same");
      cg5->Draw("same");
      std::cout<<std::endl;

    theApp.Run();
    return 0;
  }
}
