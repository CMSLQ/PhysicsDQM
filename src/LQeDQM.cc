/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/11/02 22:18:12 $
 *  $Revision: 1.1 $
 *  \author Ellie Twedt, University of Maryland 
 */

#include "Leptoquarks/PhysicsDQM/interface/LQeDQM.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Physics Objects
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
// Trigger stuff
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"

#include <vector>

#include <string>
#include <cmath>
using namespace std;
using namespace edm;
using namespace reco;



LQeDQM::LQeDQM(const ParameterSet& parameters) {
  // Get parameters from configuration file
  theElecTriggerPathToPass    = parameters.getParameter<string>("elecTriggerPathToPass");
  jetPtCut                    = parameters.getParameter<double>("jetPtCut");
  elePtCut                    = parameters.getParameter<double>("elePtCut");
  metCut                    = parameters.getParameter<double>("metCut");
  theTriggerResultsCollection = parameters.getParameter<InputTag>("triggerResultsCollection");
  theElectronCollectionLabel  = parameters.getParameter<InputTag>("electronCollection");
  theCaloJetCollectionLabel   = parameters.getParameter<InputTag>("caloJetCollection");
  theCaloMETCollectionLabel   = parameters.getParameter<InputTag>("caloMETCollection");
}

LQeDQM::~LQeDQM() { 
}


void LQeDQM::beginJob(EventSetup const& iSetup) {

  logTraceName = "LQeAnalyzer";

  LogTrace(logTraceName)<<"Parameters initialization";
  theDbe = Service<DQMStore>().operator->();
  theDbe->setCurrentFolder("Physics/ExoDQM");  // Use folder with name of PAG

  const float pi = 3.14159265;

  // Keep the number of plots and number of bins to a minimum!
  h_e1j_invMass   = theDbe->book1D("h_e1j_invMass",   "e1j Invariant Mass;InvMass (Gev)"        , 50, 0.0, 1000.0 );
  h_e2j_invMass   = theDbe->book1D("h_e2j_invMass",   "e2j Invariant Mass;InvMass (Gev)"        , 50, 0.0, 1000.0 );
  h_ee_invMass   = theDbe->book1D("h_ee_invMass",   "ee Invariant Mass;InvMass (Gev)"        , 50, 0.0, 1000.0 );
  h_eMET_transMass   = theDbe->book1D("h_eMET_transMass",   "elec+MET Transverse Mass;Transverse Mass (Gev)"        , 50, 0.0, 1000.0 );
  h_ST           = theDbe->book1D("h_ST","ST = pT 2 leading ele + pT 2 leading jets",100,0,2000);
  h_STmet           = theDbe->book1D("h_STmet","STmet = pT leading ele + pT 2 leading jets + MET;STmet (GeV)",100,0,2000);
  h_metj_transMass   = theDbe->book1D("h_metj_transMass",   "Transverse Mass MET+jet (Gev)"        , 50, 0.0, 1000.0 );
  h_jet1_et       = theDbe->book1D("h_jet1_et",       "Jet with highest E_{T} (from "+theCaloJetCollectionLabel.label()+");E_{T}(1^{st} jet) (GeV)",    20, 0., 500.0);
  h_jet2_et      = theDbe->book1D("h_jet2_et",      "Jet with 2^{nd} highest E_{T} (from "+theCaloJetCollectionLabel.label()+");E_{T}(2^{nd} jet) (GeV)",    20, 0., 500.0);
  h_jet1_eta       = theDbe->book1D("h_jet1_eta", "#eta of Leading Jet;#eta"                , 20, -5.0 , 5.0);
  h_jet2_eta       = theDbe->book1D("h_jet2_eta", "#eta of Second Jet;#eta"                 , 20, -5.0 , 5.0);
  h_jet_count    = theDbe->book1D("h_jet_count",    "Number of "+theCaloJetCollectionLabel.label()+" (E_{T} > jet Pt cut);Number of Jets", 11, -0.5, 10.5);
  h_e1_et        = theDbe->book1D("h_e1_et",  "E_{T} of Leading Electron;E_{T} (GeV)"        , 20,  0.0 , 500.0);
  h_e2_et        = theDbe->book1D("h_e2_et",  "E_{T} of Second Electron;E_{T} (GeV)"         , 20,  0.0 , 500.0);
  h_e1_eta       = theDbe->book1D("h_e1_eta", "#eta of Leading Electron;#eta"                , 20, -4.0 , 4.0);
  h_e2_eta       = theDbe->book1D("h_e2_eta", "#eta of Second Electron;#eta"                 , 20, -4.0 , 4.0);
  h_e1_phi       = theDbe->book1D("h_e1_phi", "#phi of Leading Electron;#phi"                , 22, (-1.-1./10.)*pi, (1.+1./10.)*pi );
  h_ele_count    = theDbe->book1D("h_ele_count", "Number of ele above Pt Cut", 5, -0.5, 4.5);
  h_met          = theDbe->book1D("h_met",        "Missing E_{T}; GeV"                       , 20,  0.0 , 500);
}


void LQeDQM::analyze(const Event& iEvent, const EventSetup& iSetup) {

  LogTrace(logTraceName)<<"Analysis of event # ";
  // Did it pass certain HLT path?
  Handle<TriggerResults> HLTresults;
  iEvent.getByLabel(theTriggerResultsCollection, HLTresults); 
  if ( !HLTresults.isValid() ) return;
  HLTConfigProvider hltConfig;
  hltConfig.init("HLT");
  //unsigned int triggerIndex_elec = hltConfig.triggerIndex(theElecTriggerPathToPass);
  bool passed_electron_HLT = true;
  //if (triggerIndex_elec < HLTresults->size()) passed_electron_HLT = HLTresults->accept(triggerIndex_elec);
  //if ( !(passed_electron_HLT) ) return;

  ////////////////////////////////////////////////////////////////////////////////
  // grab "gaussian sum fitting" electrons
  Handle<GsfElectronCollection> electronCollection;
  iEvent.getByLabel(theElectronCollectionLabel, electronCollection);
  if ( !electronCollection.isValid() ) return;

  // Find the highest and 2nd highest electron
  float electron_et   = -8.0;
  float electron_eta  = -8.0;
  float electron_phi  = -8.0;
  float electron2_et  = -9.0;
  float electron2_eta = -9.0;
  float electron2_phi = -9.0;
  int ele_count = 0;
  TLorentzVector e1, e2;

  // If it passed electron HLT and the collection was found, find electrons near Z mass
  if( passed_electron_HLT ) {

    for (reco::GsfElectronCollection::const_iterator recoElectron=electronCollection->begin(); recoElectron!=electronCollection->end(); recoElectron++){

      // Require electron to pass some basic cuts
      if ( recoElectron->et() < elePtCut || fabs(recoElectron->eta())>2.5 ) continue;
      ele_count++;

      if (recoElectron->et() > electron_et){
	electron2_et  = electron_et;  // 2nd highest gets values from current highest
	electron2_eta = electron_eta;
        electron2_phi = electron_phi;
        electron_et   = recoElectron->et();  // 1st highest gets values from new highest
        electron_eta  = recoElectron->eta();
        electron_phi  = recoElectron->phi();
        e1 = TLorentzVector(recoElectron->momentum().x(),recoElectron->momentum().y(),recoElectron->momentum().z(),recoElectron->p());
      } else if (recoElectron->et() > electron2_et) {
	electron2_et  = recoElectron->et();
        electron2_eta = recoElectron->eta();
        electron2_phi = recoElectron->phi();
	e2 = TLorentzVector(recoElectron->momentum().x(),recoElectron->momentum().y(),recoElectron->momentum().z(),recoElectron->p());
      }
    } // end of loop over electrons

  } // end of "are electrons valid"
  ////////////////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////////////
  // Find the highest et jet
  Handle<CaloJetCollection> caloJetCollection;
  iEvent.getByLabel (theCaloJetCollectionLabel,caloJetCollection);
  if ( !caloJetCollection.isValid() ) return;

  float jet1_et    = -8.0;
  float jet1_eta   = -8.0;
  float jet1_phi   = -8.0;
  int   jet_count = 0;
  float jet2_et   = -9.0;
  float jet2_eta  = -9.0;
  float jet2_phi  = -9.0;
  TLorentzVector j1,j2;

  for (CaloJetCollection::const_iterator i_calojet = caloJetCollection->begin(); i_calojet != caloJetCollection->end(); i_calojet++) {

    float jet_current_et = i_calojet->et();

    // if it overlaps with electron, it is not a jet
    if ( electron_et>0.0 && fabs(i_calojet->eta()-electron_eta ) < 0.2 && calcDeltaPhi(i_calojet->phi(), electron_phi ) < 0.2) continue;
    if ( electron2_et>0.0&& fabs(i_calojet->eta()-electron2_eta) < 0.2 && calcDeltaPhi(i_calojet->phi(), electron2_phi) < 0.2) continue;

    // if it has too low Et, throw away
    if (jet_current_et < jetPtCut) continue;

    jet_count++;
    if (jet_current_et > jet1_et) {
      jet2_et  = jet1_et;  // 2nd highest jet get's et from current highest
      jet2_eta = jet1_eta;
      jet2_phi = jet1_phi;
      jet1_et   = i_calojet->et(); // current highest jet gets et from the new highest
      jet1_eta  = i_calojet->eta();
      jet1_phi  = i_calojet->phi();
      j1 = TLorentzVector(i_calojet->momentum().x(),i_calojet->momentum().y(),i_calojet->momentum().z(),i_calojet->p());
    } else if (jet_current_et > jet2_et) {
      jet2_et  = i_calojet->et();
      jet2_eta = i_calojet->eta();
      jet2_phi = i_calojet->phi();
      j2 = TLorentzVector(i_calojet->momentum().x(),i_calojet->momentum().y(),i_calojet->momentum().z(),i_calojet->p());
    }
  }

  ////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  //Missing ET
  Handle<CaloMETCollection> caloMETCollection;
  iEvent.getByLabel(theCaloMETCollectionLabel, caloMETCollection);
  if ( !caloMETCollection.isValid() ) return;
  float missing_et = caloMETCollection->begin()->et();
  TLorentzVector MET_vector;
  //  MET_vector = TLorentzVector(caloMETCollection->begin()->momentum().x(),caloMETCollection->begin()->momentum().y(),caloMETCollection->begin()->momentum().z(),caloMETCollection->begin()->p());

  TVector2 v_MET, v_jet1, v_jet2, v_ele1;
  v_MET.SetMagPhi(1,caloMETCollection->begin()->phi());
  v_jet1.SetMagPhi(1,jet1_phi);
  v_jet2.SetMagPhi(1,jet2_phi);
  v_ele1.SetMagPhi(1,electron_phi);
  float deltaphi1 = v_MET.DeltaPhi(v_jet1);
  float deltaphi2 = v_MET.DeltaPhi(v_jet2);
  double MTnj1 = sqrt(2 * missing_et * jet1_et * (1 - cos(deltaphi1)) );
  double MTnj2 = sqrt(2 * missing_et * jet2_et * (1 - cos(deltaphi2)) );
  double MTne1 = sqrt(2 * electron_et * missing_et * (1 - cos(v_MET.DeltaPhi(v_ele1)) ));

  //////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  //                 Fill Histograms                                            //
  ////////////////////////////////////////////////////////////////////////////////

  bool fill_e1  = false;
  bool fill_e2  = false;

  TLorentzVector e1j1, e1j2, e2j1, e2j2, ee;
  if (electron_et > 0.0){
    fill_e1 = true;
    if (jet1_et > 0.0) {
      e1j1 = e1 + j1;
      h_e1j_invMass->Fill(e1j1.M());
    }
    if (jet2_et > 0.0) {
      e1j2 = e1 + j2;
      h_e1j_invMass->Fill(e1j2.M());
      h_STmet->Fill(jet1_et + jet2_et + electron_et + missing_et);
    }
  }
  if (electron2_et>0.0) {
    fill_e2 = true;
    ee = e1 + e2;
    h_ee_invMass->Fill(ee.M());
    if (jet1_et > 0.0) {
      e2j1 = e2 + j1;
      h_e2j_invMass->Fill(e2j1.M());
    }
    if (jet2_et > 0.0) {
      e2j2 = e2 + j2;
      h_e2j_invMass->Fill(e2j2.M());
      h_ST->Fill(jet1_et + jet2_et + electron_et + electron2_et);
    }
  }

  h_ele_count->Fill(ele_count);
  h_jet_count->Fill(jet_count);

  if (jet1_et>0.0) {
    h_jet1_et   ->Fill(jet1_et);
    h_jet1_eta  ->Fill(jet1_eta);
  }

  if (jet2_et>0.0) {
    h_jet2_et   ->Fill(jet2_et);
    h_jet2_eta   ->Fill(jet2_eta);
  }

  if (fill_e1) {
    h_e1_et      ->Fill(electron_et);
    h_e1_eta     ->Fill(electron_eta);
    h_e1_phi     ->Fill(electron_phi);
  }
  if (fill_e2) {
    h_e2_et      ->Fill(electron2_et);
    h_e2_eta     ->Fill(electron2_eta);
  }

  if (missing_et>metCut) {
    if (fill_e1)  h_eMET_transMass->Fill(MTne1);
    if (jet1_et>0.0){
      h_metj_transMass->Fill(MTnj1);
    }
    if (jet2_et>0.0) {
      h_metj_transMass->Fill(MTnj2);
    }
    h_met        ->Fill(missing_et);
  }
  ////////////////////////////////////////////////////////////////////////////////
}


void LQeDQM::endJob(void) {}

// This always returns only a positive deltaPhi
double LQeDQM::calcDeltaPhi(double phi1, double phi2) {

  double deltaPhi = phi1 - phi2;

  if (deltaPhi < 0) deltaPhi = -deltaPhi;

  if (deltaPhi > 3.1415926) {
    deltaPhi = 2 * 3.1415926 - deltaPhi;
  }

  return deltaPhi;
}
