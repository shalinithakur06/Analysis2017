
///////////////////////
// Muon Channel
///////////////////////

#include "Analyzer.h"
#include <map>

using namespace std;
void Analyzer::CutFlowAnalysis(TString url, bool isMuChannel, bool isEleChannel, TFile *outFile_){
  //check if the input file is MC or Data  
  Reader *evR_;  
  evR_ = new Reader();
  TFile *f_ = TFile::Open(url);
  int nEntries = evR_->AssignEventTreeFrom(f_);
  MyEvent *ev_;
  ev_ = evR_->GetNewEvent(1);
  TString chName = "Lepton";
  if(isMuChannel) chName = "Muon";
  if(isEleChannel) chName = "Electron";

  CutFlowProcessor(url,  chName+"/base", outFile_);
  
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  /*if(!ev_->isData){ 
    CutFlowProcessor(url,  chName+"/JESPlus",       outFile_);
    CutFlowProcessor(url,  chName+"/JESMinus",      outFile_);
    CutFlowProcessor(url,  chName+"/JERPlus",       outFile_);
    CutFlowProcessor(url,  chName+"/JERMinus",      outFile_);
    CutFlowProcessor(url,  chName+"/bTagPlus",      outFile_);
    CutFlowProcessor(url,  chName+"/bTagMinus",     outFile_);
    CutFlowProcessor(url,  chName+"/TopPtPlus",     outFile_);
    CutFlowProcessor(url,  chName+"/TopPtMinus",    outFile_);
  }
 */
  f_->Close();
  delete f_;
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void Analyzer::CutFlowProcessor(TString url, TString cutflowType, TFile *outFile_){
  int input_count_PreSel = 0;
  int input_count_ZTag = 0;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  bool isMuChannel = false;
  bool isEleChannel = false;
  if(cutflowType.Contains("Muon")) isMuChannel = true;
  if(cutflowType.Contains("Electron")) isEleChannel = true;

  //Uncertainty variations, JES, JER, MET unclustered, bTag
  int jes = 0, jer = 0, metuc = 0, bScale = 0;
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("bTagPlus"))bScale = 1;
  else if (cutflowType.Contains("bTagMinus"))bScale = -1; 
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  ////double lumiTotal = 41900;
  double lumiTotal = 41525;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  fillHisto(outFile_, cutflowType, "", "totalEvents", 10, 0, 10000000000, initialEvents, 1 );
  MyEvent *ev;
 
  //---------------------------------------------------//
  //BTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & bTagCSVfile 	= "data/CSVv2_Moriond17_B_H.csv";
  const std::string & bTagName 		= "CSVv2";
  const std::string & bTagSys 		= "central"; 
  if(bScale==1) const std::string &bTagSys 		= "up"; 
  if(bScale==-1)const std::string &bTagSys 		= "down"; 
  const std::vector<std::string> & otherSysTypes = {"up", "down"};
  //b-quark
  BTagCalibrationReader readBTagCSV_bT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readBTagCSV_cT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readBTagCSV_lT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "incl", bTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getBTagEffHistos(f);
  TString histPath("myMiniTreeProducer/MCINFO/");
  TH2D* h2_BTagEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_b"));
  TH2D* h2_BTagEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_c"));
  TH2D* h2_BTagEff_Denom_udsg 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_udsg")); 
  TH2D* h2_BTagEff_Num_bT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_bT"));
  TH2D* h2_BTagEff_Num_cT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_cT"));
  TH2D* h2_BTagEff_Num_udsgT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_udsgT")); 
 
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  double kfCount = 0;
  double n_negEvt = 0.0;
  double n_posEvt = 0.0;
  double n_oneTrig = 0.0;
  double n_twoTrig = 0.0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
    if(i > 10000) break;
    //PU histogram for MC
   /* if(i==0 && !ev->isData){
    TH1F* intimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/intimepu"));
    TH1F* trueintimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueintimepu"));
    outFile_->cd();
    intimepu->Write();
    trueintimepu->Write();
    } */
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    double genWeight = 0.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      TString sampleName_(sampleName);
      if(sampleName_.Contains("DYJetsToLL")){
        double sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight;
        evtWeight *= ev->sampleInfo.gen_weight;
	genWeight = ev->sampleInfo.gen_weight;
        if(ev->sampleInfo.gen_weight ==1) n_posEvt++;
        else n_negEvt++;
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      }	  
      //lumi weight
      else {
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight; 
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      } 
      //pileup weight
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
        float npu = pu[0];
        double weightPU = LumiWeights_.weight(npu);
        evtWeight *= weightPU;  
        fillHisto(outFile_, cutflowType, "", "puWeight", 1000, 0, 100, weightPU, 1 );
      }
      if(i==0){
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
	//cout << "Total weight of sample=" << sampleWeight << endl;
	fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents*sampleWeight);
      }
    } 
    else{ 
      if(i==0)fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents);
    }
    double topPtWt = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      if(sampleName.find("TT") != string::npos){
        vector<double>topptweights = ev->sampleInfo.topPtWeights;
        if(topptweights.size() > 0){
          topPtWt = topptweights[0]; 
          if(cutflowType.Contains("TopPtPlus")){
            topPtWt = topptweights[0];
            topPtWt = topPtWt*topPtWt;
          }
          else if(cutflowType.Contains("TopPtMinus"))
            topPtWt = 1.0;
        }
      }
    }
    fillHisto(outFile_, cutflowType, "", "SF_topPtWeights", 1000, 0, 3, topPtWt, 1 );
    evtWeight *= topPtWt; //Multiply to the total weights
        
    //---------------------------------------------------//
    //Apply triggers
    //---------------------------------------------------//
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6
    bool passTrig = false;
    bool passOneTrig = false;
    vector<string> trig = ev->hlt;
    //Muon channel
    if(isMuChannel){
      for(size_t it = 0; it < trig.size(); it++){
        if(trig[it].find("HLT_Mu50") != string::npos || 
      	      trig[it].find("HLT_TkMu50") != string::npos) {
          passTrig = true;
        }
      } 
      for(size_t it = 0; it < trig.size(); it++){
        if(trig[it].find("HLT_Mu50") != string::npos) passOneTrig = true;
      } 
      if(passTrig)n_twoTrig++;
      if(passOneTrig) n_oneTrig++;
    }
    //Electron channel
    if(isEleChannel){
      for(size_t it = 0; it < trig.size(); it++){
        //if(trig[it].find("HLT_DoubleEle33_CaloIdL") != string::npos) {
	//if(trig[it].find("HLT_DoubleEle33_CaloIdL_MW") != string::npos) {
        if(trig[it].find("HLT_Ele27_WPTight_Gsf") != string::npos) {
          passTrig = true;
        }
      }
    }
    if(!passTrig) continue;
    double nCutPass = 1.0;
    fillHisto(outFile_, cutflowType, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
   
    //---------------------------------------------------//
    //get all objets e.g. leptons, jets, vertices etc.
    //---------------------------------------------------//
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    //preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0]);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    //remove those jets which are closer to a lepton
    vector<int> j_final; j_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_init, 0.4);
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    int nLepton = 0;
    if(isMuChannel) nLepton = m_init.size();
    if(isEleChannel) nLepton = e_init.size();
    if(nLepton < 2)continue;
    double pri_vtxs = Vertices[0].totVtx;
    int charge1 = 0;
    int charge2 = 0;
    int lepton1 = 0;
    int lepton2 = 0;
    if(isMuChannel){
      lepton1 = m_init[0];
      lepton2 = m_init[1];
      charge1 = pfMuons[lepton1].charge;
      charge2 = pfMuons[lepton2].charge;
      fillHisto(outFile_, cutflowType, "", "noCharge", 10, -2, 2, 1, evtWeight);
      if(charge1 == charge2) 
        fillHisto(outFile_, cutflowType, "", "sameCharge", 10, -2, 2, 1, evtWeight);
      if(charge1 != charge2)
        fillHisto(outFile_, cutflowType, "", "oppCharge", 10, -2, 2, 1, evtWeight);
    }
    if(isEleChannel){
      lepton1 = e_init[0];
      lepton2 = e_init[1];
      charge1 = pfElectrons[lepton1].charge;
      charge2 = pfElectrons[lepton2].charge;
      fillHisto(outFile_, cutflowType, "", "noCharge", 10, -2, 2, 1, evtWeight);
      if(charge1 == charge2) 
        fillHisto(outFile_, cutflowType, "", "sameCharge", 10, -2, 2, 1, evtWeight);
      if(charge1 != charge2)
        fillHisto(outFile_, cutflowType, "", "oppCharge", 10, -2, 2, 1, evtWeight);
    }
    //both lepton should have opposite charge
    if(charge1 == charge2) continue;
    //veto third loose lepton
    bool isVeto = false;
    if(isMuChannel){
        isVeto = looseElectronVeto(-1, -1, pfElectrons) || 
            looseMuonVeto(lepton1, lepton2, pfMuons) ;
    }
    if(isEleChannel){
        isVeto = looseElectronVeto(lepton1, lepton2, pfElectrons) || 
            looseMuonVeto(-1, -1, pfMuons) ;
    }
    if(isVeto) continue;
    ///if(isMuChannel && pfMuons[lepton1].p4.pt() < 53) continue;
    if(isEleChannel && pfElectrons[lepton1].p4.pt() < 50) continue;
    //---------------------------------------------------//
    //apply lepton SF to eventWeights 
    //---------------------------------------------------//
    double leptonSF = 1.0;
    if(isMuChannel && !ev->isData){
      double lumi_BCDEF = 41.525;
      double lumi = lumi_BCDEF;
      //get muon scale factor for fist muon
      //trigger 	
      double muSFtrig_BCDEF1 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFtrig1 		= (muSFtrig_BCDEF1); 
      //identification
      double muSFid_BCDEF1 	= getMuonSF(h2_idSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFid1 		= (muSFid_BCDEF1); 
      //isolation 
      double muSFiso_BCDEF1 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFiso1 		= (muSFiso_BCDEF1); 
      double muSF1 = muSFtrig1*muSFid1*muSFiso1;	
      
      //get muon scale factor for 2nd muon
      //trigger 	
      double muSFtrig_BCDEF2 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFtrig2 		= (muSFtrig_BCDEF2); 
      //identification
      double muSFid_BCDEF2 	= getMuonSF(h2_idSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFid2 		= (muSFid_BCDEF2); 
      //isolation 
      double muSFiso_BCDEF2 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFiso2 		= (muSFiso_BCDEF2); 
      //tracking 
      double muSF2 = muSFtrig2*muSFid2*muSFiso2;	
      leptonSF = muSF1*muSF2;
      ///cout << "leptonSF = " << muSF1*muSF2 << endl;
    } 
    if(isEleChannel && !ev->isData){
      ///double eleSF1 =0;
      double eleSF1=1.0; double ele_heepEB_SF1=1.0; double ele_heepEE_SF1=1.0;
      double ele_trigSF1=1.0;
      ele_trigSF1=0.991; ele_heepEB_SF1=0.967; ele_heepEE_SF1=0.973;
      float e1Eta = TMath::Abs(pfElectrons[lepton1].p4.eta());
      bool isEB = e1Eta<=1.444;
      ///bool isEE = e1Eta>1.566;
      ///bool isEB = false;
      ///if(e1Eta<=1.444) isEB=true;
      double ele_recoSF1       = getEleSF(h2_ele_recoSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      if(isEB){
	eleSF1 = ele_recoSF1*ele_trigSF1*ele_heepEB_SF1;
	cout<<"isEB="<< isEB <<endl;
	cout<<"ele SF1 = "<< eleSF1 <<endl;
      }	
      else {
	   eleSF1 = ele_recoSF1*ele_trigSF1*ele_heepEE_SF1;
	   cout<<"isEB="<< isEB <<endl;
	   cout<<"ele SF1 ="<< eleSF1 <<endl;
      }  
      //SFs for 2nd electron
      double eleSF2 =1.0; double ele_heepEB_SF2=1.0; double ele_heepEE_SF2=1.0;
      double ele_trigSF2=1.0;
      ele_trigSF2=0.991; ele_heepEB_SF2=0.967; ele_heepEE_SF2=0.973;
      float e2Eta = TMath::Abs(pfElectrons[lepton2].p4.eta());
      double ele_recoSF2       = getEleSF(h2_ele_recoSF, pfElectrons[lepton2].eleSCEta, pfElectrons[lepton2].p4.pt());
      if(isEB){
       eleSF2 = ele_recoSF2*ele_trigSF2*ele_heepEB_SF2;
       cout<<"isEB="<< isEB <<endl;
       cout<<"ele SF2 ="<< eleSF2 << endl;
       }
       else {
         eleSF2 = ele_recoSF2*ele_trigSF2*ele_heepEE_SF2;
         //cout<<"isEB="<< isEB <<endl;
         //cout<<"ele SF2 ="<< eleSF2 << endl;
       }
       leptonSF = eleSF1*eleSF2;	       
       cout<<"Lepton SF is ="<< leptonSF << endl;
       //combined SF (Scale factors are applied on MC only.)
       /*double eleSF =1.0;
       if(!ev->isData){
         eleSF = eleSF1*eleSF2;
         leptonSF = eleSF1*eleSF2;	       
         ////cout << "leptonSF = " << eleSF2 << endl;      
         ////cout << "leptonSF = " << eleSF1*eleSF2 << endl;      
       } */

      //This is cut-based ID, we are using Heep ID
      //double ele_medium_idSF1  = getEleSF(h2_ele_medium_idSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      ///double ele_trigSF1       = getEleTrigSF(h2_ele_trigSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      ////double ele_heep_SF1      = getEleHeep2SF(tg_heep_SF, pfElectrons[lepton1].eleSCEta);
      ///eleSF1 = ele_recoSF1*ele_trigSF1*ele_heep_SF1;  

      ///double eleSF2 =0;
      ///double ele_recoSF2       = getEleSF(h2_ele_recoSF, pfElectrons[lepton2].eleSCEta, pfElectrons[lepton2].p4.pt());
      //double ele_medium_idSF2  = getEleSF(h2_ele_medium_idSF, pfElectrons[lepton2].eleSCEta, pfElectrons[lepton2].p4.pt());
      ///double ele_trigSF2       = getEleTrigSF(h2_ele_trigSF, pfElectrons[lepton2].eleSCEta,  pfElectrons[lepton2].p4.pt());
      ///double ele_heep_SF2      = getEleHeep2SF(tg_heep_SF, pfElectrons[lepton2].eleSCEta);
      ///eleSF2 = ele_recoSF2*ele_trigSF2*ele_heep_SF2;
      ///double eleSF = 1.0;
      ///leptonSF = eleSF1*eleSF2;
    }
    evtWeight *= leptonSF;
    cout << "leptonSF = " << evtWeight << endl;      
    
    fillHisto(outFile_, cutflowType, "", "leptonSF", 1000, 0, 100, leptonSF, 1 );
    string cutflowType_(cutflowType);
    cutflowType_ = cutflowType;
    
    //---------------------------------------------------//
    //get 4 vector for Z boson
    //---------------------------------------------------//
    MyLorentzVector vZ; 
    if(isMuChannel) vZ = pfMuons[lepton1].p4 + pfMuons[lepton2].p4;
    if(isEleChannel) vZ = pfElectrons[lepton1].p4 + pfElectrons[lepton2].p4;
    if(j_final.size()==0) continue;
    int count_jets = j_final.size();
    bool isControlSel = true;
    if(vZ.mass() < 60) isControlSel = false;
    //---------------------------------------------------//
    //Fill histos with for Control Plots
    //---------------------------------------------------//
    //fill histos for muon
    double leptonPt1 = 0.0;
    double leptonPt2 = 0.0; 
    if(isMuChannel){
      leptonPt1 = pfMuons[lepton1].p4.pt();
      leptonPt2 = pfMuons[lepton2].p4.pt(); 
    }
    if(isEleChannel){
      leptonPt1 = pfElectrons[lepton1].p4.pt();
      leptonPt2 = pfElectrons[lepton2].p4.pt(); 
    }
    double metPt = 0; 
    metPt = metWithJESJER(pfJets, &j_final, met, jes, jer, ev->isData);
    if(isControlSel){
      double dR1 = 0.0;
      double dR2 = 0.0;
      double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
      double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
      //fill histos for jets
      double count_bJetT = 0.0;
      ///double count_CSVT_SF = 0.0;
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
        if(isMuChannel){
          dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton1].p4);
          dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton2].p4);
        }
        if(isEleChannel){
          dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton1].p4);
          dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton2].p4);
        }
        fillHisto(outFile_, cutflowType_, "ControlP","dR1", 50, 0, 10, dR1, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","dR2", 50, 0, 10, dR2, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","dR", 50, 0, 10, dR1, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","dR", 50, 0, 10, dR2, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","pt_jet", 250, 0, 10000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","mass_jet", 250, 0, 5000, pfJets[ind_jet].p4.M(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","ak8Pmass", 250, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
        pfCISV = pfJets[ind_jet].bDiscriminator["pfDeepCSVb"];
        pfCMVA = pfJets[ind_jet].bDiscriminator["pfDeepCSVbb"];
        fillHisto(outFile_, cutflowType_, "ControlP", "pfCISV", 100, -2, 2, pfCISV, evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
        if(pfCISV + pfCMVA  > 0.8001){  //0.8484//tight working point
          count_bJetT++;
	  ///count_CSVT_SF++; 
          fillHisto(outFile_, cutflowType_, "ControlP", "pt_bjet", 100, 0, 1000, jetPt, evtWeight );
          fillHisto(outFile_, cutflowType_, "ControlP", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
        }
      }
      fillHisto(outFile_, cutflowType_, "ControlP","multi_bjet",  16, 0.5, 15.5, count_bJetT, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","MET",  250, 0, 10000, metPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","multiLep",  15, 0.5, 15.5, nLepton, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","pt1stLep", 250, 0, 10000, leptonPt1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","pt2ndLep", 250, 0, 10000, leptonPt2, evtWeight );
      fillHisto2D(outFile_, cutflowType_,"ControlP", "ptLep1_ptLep2", 250, 0, 10000, leptonPt1,250, 0, 10000, leptonPt2, 1);
      double etaLep1 = 0.0;
      double etaLep2 =0.0;
      if(isMuChannel){
        etaLep1 =  pfMuons[lepton1].p4.eta();
        etaLep2 =  pfMuons[lepton2].p4.eta();
      }
      if(isEleChannel){
        etaLep1 =  pfElectrons[lepton1].p4.eta();
        etaLep2 =  pfElectrons[lepton2].p4.eta();
      }
      fillHisto(outFile_, cutflowType_, "ControlP","eta1stLep", 50, -5, 5, etaLep1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","eta2ndLep", 50, -5, 5, etaLep2, evtWeight );
      
      //fill histos for Z boson
      fillHisto(outFile_, cutflowType_, "ControlP","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","mll", 250, 0, 10000, vZ.M(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","pt_Z",  250, 0, 10000, vZ.Pt(), evtWeight );
      /*if(vZ.M() > 80 && vZ.M()< 100){
        fillHisto(outFile_, cutflowType_, "ControlP2","pt_Z",  200, 0, 1000, vZ.Pt(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP2","pz_Z",  200, 0, 1000, vZ.Pz(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP2","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
        fillHisto(outFile_, cutflowType_, "ControlP2","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
      }
      fillHisto(outFile_, cutflowType_, "ControlP2","mll", 20, 50, 150, vZ.M(), evtWeight ); */
      //fill histos for nvtx
      fillHisto(outFile_, cutflowType_, "ControlP","nvtx", 50, 0, 100, pri_vtxs, evtWeight );
      for(std::size_t n=0; n<Vertices.size(); n++){
        fillHisto(outFile_, cutflowType_, "ControlP","rhoAll", 50, 0, 100, Vertices[n].rhoAll, evtWeight );
      }
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    }

    //---------------------------------------------------//
    //Fill histos with pre-selection
    // Impact of PreSel cut: 
    // Choose only those events where ALL the jets satisfy 
    // the below jet selection criteria
    //---------------------------------------------------//
    std::vector<int> storePreSel;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      double jetEta = fabs(pfJets[ind_jet].p4.eta());
      double jetPmass = pfJets[ind_jet].ak8Pmass;
      double dR1 = 0.0;
      double dR2 = 0.0;
      if(isMuChannel){
        dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton1].p4);
        dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton2].p4);
      }
      if(isEleChannel){
        dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton1].p4);
        dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton2].p4);
      }
      if(jetPt > 100 && jetEta < 2.5 && dR1 > 0.8 && 
		      dR2 > 0.8 && vZ.Pt() > 100 && jetPmass > 40)storePreSel.push_back(1);
      else storePreSel.push_back(0);
    }
    //check pre selection condition 
    bool isPreSel = true;
    for(unsigned int f = 0; f <storePreSel.size(); f++){
      if (storePreSel[f]==0) isPreSel=false;
    }
    if(isControlSel && isPreSel){
      double pfCISV = 0.0; 
      double pfCMVA = 0.0;
      double count_bJetT = 0.0;
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
        double jetEta = fabs(pfJets[ind_jet].p4.eta());
        double jetPmass = pfJets[ind_jet].ak8Pmass;
        double dR1 = 0.0;
        double dR2 = 0.0;
        if(isMuChannel){
          dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton1].p4);
          dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[lepton2].p4);
        }
        if(isEleChannel){
          dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton1].p4);
          dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[lepton2].p4);
        }
        fillHisto(outFile_, cutflowType_, "PreSel","dR1", 50, 0, 10, dR1, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","dR2", 50, 0, 10, dR2, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","dR", 50, 0, 10, dR1, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","dR", 50, 0, 10, dR2, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","pt_jet", 250, 0, 10000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","mass_jet", 250, 0, 5000, pfJets[ind_jet].p4.M(), evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","ak8Pmass", 250, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
	pfCISV = pfJets[ind_jet].bDiscriminator["pfDeepCSVb"];
	pfCMVA = pfJets[ind_jet].bDiscriminator["pfDeepCSVbb"];
        fillHisto(outFile_, cutflowType_, "PreSel", "pfCISV", 100, -2, 2, pfCISV, evtWeight );
        fillHisto(outFile_, cutflowType_, "PreSel", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
	if(pfCISV +pfCMVA > 0.8001){
	  count_bJetT++;
	  fillHisto(outFile_, cutflowType_, "PreSel", "pt_bjet", 100, 0, 1000, jetPt, evtWeight );
	  fillHisto(outFile_, cutflowType_, "PreSel", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
       }	  
     }
      fillHisto(outFile_, cutflowType_, "PreSel","multi_bjet", 16, -0.5, 15.5, count_bJetT, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","MET",  250, 0, 10000, metPt, evtWeight );
      //fillHisto(outFile_, cutflowType_, "PreSel","pfJets_size", 15, 0, 15, pfJets.size(), evtWeight );
      ///if(muonPt1 <100) continue;    
      //fill histos for muon
      fillHisto(outFile_, cutflowType_, "PreSel","multiLep",  15, 0.5, 15.5, nLepton, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","pt1stLep", 250, 0, 10000, leptonPt1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","pt2ndLep", 250, 0, 10000, leptonPt2, evtWeight );
      fillHisto2D(outFile_, cutflowType_,"PreSel", "ptLep1_ptLep2", 250, 0, 10000, leptonPt1,250, 0, 10000, leptonPt2, 1);
      double etaLep1 = 0.0;
      double etaLep2 =0.0;
      if(isMuChannel){
        etaLep1 =  pfMuons[lepton1].p4.eta();
        etaLep2 =  pfMuons[lepton2].p4.eta();
      }
      if(isEleChannel){
        etaLep1 =  pfElectrons[lepton1].p4.eta();
        etaLep2 =  pfElectrons[lepton2].p4.eta();
      }
      fillHisto(outFile_, cutflowType_, "PreSel","eta1stLep", 50, -5, 5, etaLep1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","eta2ndLep", 50, -5, 5, etaLep2, evtWeight );
      //fill histos for Z boson
      fillHisto(outFile_, cutflowType_, "PreSel","pt_Z",  250, 0, 10000, vZ.Pt(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","mll", 250, 0, 10000, vZ.M(), evtWeight );
      //fill histos for nvtx
      fillHisto(outFile_, cutflowType_, "PreSel","nvtx", 50, 0, 100, pri_vtxs, evtWeight );
      for(std::size_t n=0; n<Vertices.size(); n++){
        fillHisto(outFile_, cutflowType_, "PreSel","rhoAll", 50, 0, 100, Vertices[n].rhoAll, evtWeight );
      }
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight);
      input_count_PreSel++;
      if(input_count_PreSel%1000==0)
      cout << "input count after PreSel: "<< input_count_PreSel << endl;
    }
    
    //---------------------------------------------------//
    // fill histo after ZTag:
    //---------------------------------------------------//
    // tag a jet as Z-boson
    vector<size_t> allZjet;
    vector<size_t> allZjetBkgEst;
    bool isZTagSel = false;
    int countFatJet = 0;
    int countFatJetBkgEst = 0;
    double count_bJetT = 0.0;
    double count_bJetTBkgEst = 0.0;
    double count_bJetM = 0.0;
    int cut1 = 0.0;
    int cut2 = 0.0;
    int cut3 = 0.0;
    if(isControlSel && isPreSel){
    double pfCISV = 0.0;
    double pfCMVA = 0.0;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      double jetPmass = pfJets[ind_jet].ak8Pmass;
      double ak8Tau21 = pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1;
      pfCISV = pfJets[ind_jet].bDiscriminator["pfDeepCSVb"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfDeepCSVbb"];
      if(jetPt >200 && jetPmass > 65 && jetPmass < 95 && vZ.M()> 200 && ak8Tau21 < 0.60){
        countFatJet ++;
	allZjet.push_back(ijet);
      }
      else{ //to veto non-fat jet
	fillHisto(outFile_, cutflowType_, "ZTag", "pfCISV", 100, -2, 2, pfCISV, evtWeight );
	fillHisto(outFile_, cutflowType_, "ZTag", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
	if(pfCISV + pfCMVA > 0.8001){
          count_bJetT++;
	  fillHisto(outFile_, cutflowType_, "ZTag", "pt_bjet", 100, 0, 1000, jetPt, evtWeight );
	  fillHisto(outFile_, cutflowType_, "ZTag", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
        }
	if(pfCISV + pfCMVA > 0.8001){
	    count_bJetM++;
        }	
      }//else
      ///for bkg estimation
      if(jetPt > 200 && jetPmass > 65 && jetPmass < 95 ){
	countFatJetBkgEst ++;
        allZjetBkgEst.push_back(ijet);
      }	
      else{//to veto non-fat jet
	if(pfCISV + pfCMVA > 0.8001){
	  count_bJetTBkgEst++;
        }
      } //else
    }      
    double new_evtWeight = 1.0;
    if(countFatJet==1){
      //apply tau21 scale factor
      //https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging#2016_scale_factors_and_correctio
      double new_evtWeight = 1.0;
      if(!ev->isData) new_evtWeight = 1.136*0.97*evtWeight; 
      MyLorentzVector vZmax;
      MyLorentzVector vZmin;
      if(isMuChannel){
        vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[lepton1].p4;
        vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[lepton2].p4;
      }
      if(isEleChannel){
        vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[lepton1].p4;
        vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[lepton2].p4;
      }
      fillHisto(outFile_, cutflowType_, "ZTag","multi_bjet", 16, -0.5, 15.5, count_bJetT, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","multi_bjetM", 16, -0.5, 15.5, count_bJetM, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","MET",  250, 0, 10000, metPt, new_evtWeight );
      fillHisto2D(outFile_, cutflowType_,"ZTag", "Nbjet_MET",10, -0.5, 9.5, count_bJetT, 1000, 0, 10000, metPt, new_evtWeight );
      fillHisto2D(outFile_, cutflowType_,"ZTag", "NbjetM_MET",10, -0.5, 9.5, count_bJetM, 1000, 0, 10000, metPt, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mass_fatjet", 250, 0, 5000, pfJets[j_final[allZjet[0]]].p4.M(), new_evtWeight );      
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, new_evtWeight );
      //fill histos for muon
      fillHisto(outFile_, cutflowType_, "ZTag","multiLep",  15, 0.5, 15.5, nLepton, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","pt1stLep", 500, 0, 10000, leptonPt1, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","pt2ndLep", 500, 0, 10000, leptonPt2, new_evtWeight );
      fillHisto2D(outFile_, cutflowType_,"ZTag", "ptLep1_ptLep2", 500, 0, 10000, leptonPt1,500, 0, 10000, leptonPt2, 1);
      double etaLep1 = 0.0;
      double etaLep2 =0.0;
      if(isMuChannel){
        etaLep1 =  pfMuons[lepton1].p4.eta();
        etaLep2 =  pfMuons[lepton2].p4.eta();
      }
      if(isEleChannel){
        etaLep1 =  pfElectrons[lepton1].p4.eta();
        etaLep2 =  pfElectrons[lepton2].p4.eta();
      }
      fillHisto(outFile_, cutflowType_, "ZTag","eta1stLep", 50, -5, 5, etaLep1, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","eta2ndLep", 50, -5, 5, etaLep2, new_evtWeight );
      //fill histos for Z boson
      fillHisto(outFile_, cutflowType_, "ZTag","pt_Z",  250, 0, 10000, vZ.Pt(), new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","eta_Z",  50, -5, 5, vZ.Rapidity(), new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","phi_Z",  50, -5, 5, vZ.Phi(), new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mll", 250, 0, 10000, vZ.M(), new_evtWeight );
      //Section 5.1 of https://drive.google.com/drive/folders/1e8PiEjw7sWXJn7-ET8TaK81QOgFHV0_g
      double mlZmax = 0.0;
      double mlZmin = 0.0;
      if(vZmax.M() >= vZmin.M()){
        mlZmax = vZmax.M();
        mlZmin = vZmin.M();
      }
      else{
        mlZmax = vZmin.M();
        mlZmin = vZmax.M();
      }
      fillHisto(outFile_, cutflowType_, "ZTag","genWeight", 500, -2, 2, genWeight, 1);
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min", 250, 0, 10000, mlZmin, new_evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max", 250, 0, 10000, mlZmax, new_evtWeight );
      fillHisto2D(outFile_, cutflowType_,"ZTag", "mlZmin_mlZmax",250, 0, 10000, mlZmin, 250, 0, 10000, mlZmax, new_evtWeight);
      //fill histos for jets
      ///for(size_t ijet = 0; ijet < j_init.size(); ijet++){
        ///int ind_jet = j_init[ijet];
        double jetPt = jetPtWithJESJER(pfJets[j_final[allZjet[0]]], jes, jer);
        double dR1 = 0.0;
        double dR2 = 0.0;
        if(isMuChannel){
          dR1 = DeltaR(pfJets[j_final[allZjet[0]]].p4, pfMuons[lepton1].p4);
          dR2 = DeltaR(pfJets[j_final[allZjet[0]]].p4, pfMuons[lepton2].p4);
        }
        if(isEleChannel){
          dR1 = DeltaR(pfJets[j_final[allZjet[0]]].p4, pfElectrons[lepton1].p4);
          dR2 = DeltaR(pfJets[j_final[allZjet[0]]].p4, pfElectrons[lepton2].p4);
        }
        fillHisto(outFile_, cutflowType_, "ZTag","dR1", 50, 0, 10, dR1, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","dR2", 50, 0, 10, dR2, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","dR", 50, 0, 10, dR1, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","dR", 50, 0, 10, dR2, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","pt_jet", 250, 0, 10000, jetPt, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","eta_jet", 50, -5, 5, pfJets[j_final[allZjet[0]]].p4.eta(), new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","mass_jet", 250, 0, 5000, pfJets[j_final[allZjet[0]]].p4.M(), new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","phi_jet", 50, -5, 5, pfJets[j_final[allZjet[0]]].p4.phi(), new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","ak8Pmass", 250, 0, 5000, pfJets[j_final[allZjet[0]]].ak8Pmass, new_evtWeight );
        fillHisto(outFile_, cutflowType_, "ZTag","ak8Tau21", 50, 0, 5, pfJets[j_final[allZjet[0]]].ak8Tau2/pfJets[j_final[allZjet[0]]].ak8Tau1, new_evtWeight );
	fillHisto(outFile_, cutflowType_, "ZTag","final_multi_jet", 15, 0, 15, count_jets, new_evtWeight );
      ///}
      //fill histos for nvtx
      fillHisto(outFile_, cutflowType_, "ZTag","nvtx", 50, 0, 100, pri_vtxs, new_evtWeight );
      for(std::size_t n=0; n<Vertices.size(); n++){
        fillHisto(outFile_, cutflowType_, "ZTag","rhoAll", 50, 0, 100, Vertices[n].rhoAll, new_evtWeight );
      }
      input_count_ZTag++;
      if(input_count_ZTag%1000==0)
      cout << "input count after ZTag: "<< input_count_ZTag << endl;
    }

    //---------------------------------------------------//
    //apply b-veto and MET cut
    //---------------------------------------------------//
    if(countFatJet==1 && metPt < 200){
      nCutPass++;
      fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, new_evtWeight );
    }
    //---------------------------------------------------//
    //apply B-tagging
    //---------------------------------------------------//
    bool isBTagVeto = true; 
    if(count_bJetT >=2) isBTagVeto=false;
    //Apply b-tag SF
    //https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
    double pmc_btag = 1.0;
    double pdata_btag = 1.0;
    double bTagWt = 1.0; 
    if(countFatJet==1 && metPt < 200 && count_bJetT==0){
      if(!ev->isData){
        for(size_t ijet = 0; ijet < j_final.size(); ijet++){
           int ind_jet = j_final[ijet];
           double pMC_ = 1.0;
           double pData_ = 1.0;
           //b-quark
           if(abs(pfJets[ind_jet].partonFlavour) ==5){
             pMC_ = getBTagPmcSys(h2_BTagEff_Num_bT, h2_BTagEff_Denom_b, pfJets[ind_jet]); 
             pData_ = getBTagPdataSys(readBTagCSV_bT, h2_BTagEff_Num_bT, h2_BTagEff_Denom_b, pfJets[ind_jet],bScale);
           }
           //c-quark
           else if(abs(pfJets[ind_jet].partonFlavour) ==4){ 
             pMC_ = getBTagPmcSys(h2_BTagEff_Num_cT, h2_BTagEff_Denom_c, pfJets[ind_jet]); 
             pData_ = getBTagPdataSys(readBTagCSV_cT, h2_BTagEff_Num_cT, h2_BTagEff_Denom_c, pfJets[ind_jet],bScale);
           }
           //other quarks and gluon
           else{ 
             pMC_ = getBTagPmcSys(h2_BTagEff_Num_udsgT, h2_BTagEff_Denom_udsg, pfJets[ind_jet]); 
             pData_ = getBTagPdataSys(readBTagCSV_lT, h2_BTagEff_Num_udsgT, h2_BTagEff_Denom_udsg, pfJets[ind_jet], bScale); 
           }
           pmc_btag = pmc_btag*pMC_;
           pdata_btag = pdata_btag*pData_;
         }
       bTagWt = pdata_btag/pmc_btag;
       }
       new_evtWeight *= bTagWt;
       nCutPass++;
       fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, new_evtWeight);
       fillHisto(outFile_, cutflowType, "", "bTagWeight", 100, 0, 2, bTagWt, 1);
       nCutPass++;
       fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, new_evtWeight );
     } 

    //---------------------------------------------------//
    //for background estimation
    //---------------------------------------------------//
    bool isBkgEst = true;
    //if(countFatJetBkgEst==1){
    if(countFatJetBkgEst==1 && metPt < 200 && count_bJetTBkgEst==0){
      if(!ev->isData) evtWeight *= 1.01; 
      MyLorentzVector vZmax;
      MyLorentzVector vZmin;  
      if(isMuChannel){
        vZmax = pfJets[j_final[allZjetBkgEst[0]]].p4 + pfMuons[lepton1].p4;
        vZmin = pfJets[j_final[allZjetBkgEst[0]]].p4 + pfMuons[lepton2].p4;
      }
      if(isEleChannel){
        vZmax = pfJets[j_final[allZjetBkgEst[0]]].p4 + pfElectrons[lepton1].p4;
        vZmin = pfJets[j_final[allZjetBkgEst[0]]].p4 + pfElectrons[lepton2].p4;
      }
      double mlZmax = 0.0;
      double mlZmin = 0.0;
      if(vZmax.M() >= vZmin.M()){
        mlZmax = vZmax.M();
        mlZmin = vZmin.M();
      }
      else{
        mlZmax = vZmin.M();
        mlZmin = vZmax.M();
      }
      double ml1l2 = vZ.M();
      double tau21 = pfJets[j_final[allZjetBkgEst[0]]].ak8Tau2/pfJets[j_final[allZjetBkgEst[0]]].ak8Tau1;
      fillHisto2D(outFile_, cutflowType_,"ZTag", "ml1l2_tau21",10, 0, 1, tau21, 1000, 0, 10000, ml1l2, evtWeight);
      bool isRegionA = false;
      bool isRegionB = false;
      bool isRegionC = false;
      bool isRegionD = false;
      if(ml1l2 > 200 && tau21 < 0.60) isRegionA = true;
      if(ml1l2 < 200 && tau21 < 0.60) isRegionB = true;
      if(ml1l2 > 200 && tau21 > 0.60) isRegionC = true;
      if(ml1l2 < 200 && tau21 > 0.60) isRegionD = true;
      if(isRegionA){ //for sanity checks
        fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_ZTag", 250, 0, 10000, mlZmin, evtWeight );
      }
      //L-cut
      vector<string> sigMass;
      vector<double> lCutMin;
      vector<double> lCutMax;
      sigMass.push_back("250");   lCutMax.push_back(440);  lCutMin.push_back(300);
      sigMass.push_back("500");   lCutMax.push_back(450);  lCutMin.push_back(560);
      sigMass.push_back("750");   lCutMax.push_back(700);  lCutMin.push_back(900);
      sigMass.push_back("1000");  lCutMax.push_back(950);  lCutMin.push_back(1080);
      sigMass.push_back("1250");  lCutMax.push_back(950); lCutMin.push_back(1370);
      sigMass.push_back("1500");  lCutMax.push_back(950); lCutMin.push_back(1700);
      sigMass.push_back("1750");  lCutMax.push_back(950); lCutMin.push_back(1950);
      sigMass.push_back("2000");  lCutMax.push_back(950); lCutMin.push_back(2200);
      sigMass.push_back("2500");  lCutMax.push_back(950); lCutMin.push_back(2700);
      sigMass.push_back("3000");  lCutMax.push_back(950); lCutMin.push_back(3200);
      sigMass.push_back("3500");  lCutMax.push_back(950); lCutMin.push_back(3700);
      sigMass.push_back("4000");  lCutMax.push_back(950); lCutMin.push_back(4200);
      sigMass.push_back("4500");  lCutMax.push_back(950); lCutMin.push_back(4700);
      sigMass.push_back("5000");  lCutMax.push_back(950); lCutMin.push_back(5200);

      vector<double>yStepVec;     vector<double>xStepVec;   vector<string>dirVec;         
      yStepVec.push_back(0);      xStepVec.push_back(0);    dirVec.push_back("LCut1");
      yStepVec.push_back(50);     xStepVec.push_back(0);    dirVec.push_back("LCut2");
      yStepVec.push_back(100);    xStepVec.push_back(0);    dirVec.push_back("LCut3");
      yStepVec.push_back(200);    xStepVec.push_back(0);    dirVec.push_back("LCut4");
      yStepVec.push_back(300);    xStepVec.push_back(0);    dirVec.push_back("LCut5");
      yStepVec.push_back(400);    xStepVec.push_back(0);    dirVec.push_back("LCut6");
      yStepVec.push_back(500);    xStepVec.push_back(0);    dirVec.push_back("LCut7");
      yStepVec.push_back(600);    xStepVec.push_back(0);    dirVec.push_back("LCut8");
      yStepVec.push_back(700);    xStepVec.push_back(0);    dirVec.push_back("LCut9");
      yStepVec.push_back(800);    xStepVec.push_back(0);    dirVec.push_back("LCut10");

      for(unsigned int l =0; l<sigMass.size(); l++){
        double max = lCutMax[l];
        double min = lCutMin[l];
        TString mass = sigMass[l];
        for(unsigned int d = 0; d<dirVec.size(); d++){
          TString dir = dirVec[d];
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d])){
          fillHisto(outFile_, cutflowType_, dir+"/ABCD","ABCD_mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/ABCD","ABCD_mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/ABCD","ABCD_genWeight", 500, -2, 2, genWeight, 1);
          }
          //Region-A
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d]) && isRegionA){
          fillHisto(outFile_, cutflowType_, dir+"/A","A_mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/A","A_mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/A","A_genWeight", 500, -2, 2, genWeight, 1);
          }
          //Region-B
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d]) && isRegionB){
          fillHisto(outFile_, cutflowType_, dir+"/B","B_mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/B","B_mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/B","B_genWeight", 500, -2, 2, genWeight, 1);
          }
          //Region-C
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d]) && isRegionC){
          fillHisto(outFile_, cutflowType_, dir+"/C","C_mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/C","C_mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/C","C_genWeight", 500, -2, 2, genWeight, 1);
          }
          //Region-D
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d]) && isRegionD){
          fillHisto(outFile_, cutflowType_, dir+"/D","D_mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/D","D_mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          fillHisto(outFile_, cutflowType_, dir+"/D","D_genWeight", 500, -2, 2, genWeight, 1);
          }
          //Region-A; needed for old codes
          if(mlZmax > (max + yStepVec[d]) && mlZmin < (min + xStepVec[d]) && isRegionA){
          fillHisto(outFile_, cutflowType_, dir,"mlZ_min_sig"+mass, 250, 0, 10000, mlZmin, evtWeight );
          fillHisto(outFile_, cutflowType_, dir,"mlZ_max_sig"+mass, 250, 0, 10000, mlZmax, evtWeight );
          }
        }
      }
    }//bkg estimation loop
    //if(i > 200000) break;
  }//event loop
  cout<<"Total events  = "<<nEntries<<endl;
  cout<<"Total events with negative weight = "<<n_negEvt<<endl;
  cout<<"Total events with positive weight = "<<n_posEvt<<endl;

  double effective_evt = (n_posEvt-n_negEvt);
  double amcnlo_weight = 1.0;
  if(effective_evt !=0) amcnlo_weight = effective_evt/nEntries;
  fillHisto(outFile_, cutflowType, "", "positive_weight", 10, -2, 2, -1, n_posEvt);
  fillHisto(outFile_, cutflowType, "", "negative_weight", 10, -2, 2,  1, n_negEvt);
  fillHisto(outFile_, cutflowType, "", "amcnlo_weight", 10, 0, 1, amcnlo_weight, 1 );
  fillHisto(outFile_, cutflowType, "", "oneTrig", 10, -2, 4, 1, n_oneTrig);
  fillHisto(outFile_, cutflowType, "", "twoTrig", 10, -2, 4, 2, n_twoTrig);
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+"outputFile_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel(9);
  
  //Local 
  ///TString pathLocal = "outFile_.root";
  ///CutFlowAnalysis(pathLocal, true, false, outFile_); 
  ///CutFlowAnalysis(pathLocal, false, true, outFile_); 
  
  //T2
  
  TString pathT2 = "/cms/store/user/sthakur/ntuple_for2017MC_20191101/MC_20191101/DYJetsToLL_M50_MC_20191101/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M50_MC_20191101/191101_121023/0000/DYJetsToLL_M50_MC_20191101_Ntuple_1.root";
  ///CutFlowAnalysis("root://se01.indiacms.res.in:1094/"+pathT2, true, false, outFile_);
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/"+pathT2, false, true, outFile_);
  
  //================
  //condor submission
  //================
  ///CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", true, false, outFile_);
  ///CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", false, true, outFile_);
  
  outFile_->Write(); 
  outFile_->Close();
} 
