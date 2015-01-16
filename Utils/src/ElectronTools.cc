// $Id: ElectronTools.cc,v 1.55 2013/04/05 13:13:31 arapyan Exp $

#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"

ClassImp(mithep::ElectronTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ElectronTools::ElectronTools()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassCustomID(const Electron *ele, EElIdType idType) {

  Double_t  fCuts[6][8];             //!custom id cuts

  Double_t tightcuts[6][8]={
    {0.086, 0.1, 0.052, 0.0, 0.050, 0.059, 0.061, 0.0},   //hovere
    {0.011, 0.011, 0.011, 0.0, 0.033, 0.029, 0.030, 0.0},     //sigmaetaeta
    {0.038, 0.024, 0.045, 0.0, 0.034, 0.017, 0.026, 0.0},   //deltaphiin
    {0.0081, 0.0029, 0.0051, 0.0, 0.0070, 0.0062, 0.0088, 0.0}, //deltaetain
    {0.0,    0.9,    0.0,    0.0, 0.0,    0.78,   0.0,    0.0},             //eoverp
    {0.8,0.2,0.9,0,0,0,0,0}};                              //extra cuts fbrem and E_Over_P 

  Double_t loosecuts[6][8]={
    {0.12,  0.12,  0.12,  0.12,  0.10,   0.10,   0.10,   0.10 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.03,   0.03,   0.03,   0.03 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.007,  0.007,  0.007,  0.007}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,	 0.0,	 0.0,	 0.0  }, //eoverp
    {0.0,   0.0,   0,	  0,	 0,	 0,	 0,	 0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPointFakeable[6][8] = {
    {0.12,  0.12,  0.12,  0.12,  0.10,   0.10,   0.10,   0.10  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.15,  0.15,  0.15,  0.15,  0.10,   0.10,   0.10,   0.10  }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.009,  0.009,  0.009,  0.009 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t Veto[6][8] = {
    {0.15,  0.15,  0.15,  0.15,  999.,   999.,   999.,   999.  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.010,  0.010,  0.010,  0.010 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint95[6][8] = {
    {0.15,  0.15,  0.15,  0.15,  0.07,   0.07,   0.07,   0.07  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.010,  0.010,  0.010,  0.010 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint90[6][8] = {
    {0.12,  0.12,  0.12,  0.12,  0.05,   0.05,   0.05,   0.05  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.009,  0.009,  0.009,  0.009 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint85[6][8] = {
    {0.04,  0.04,  0.04,  0.04,  0.025,  0.025,  0.025,  0.025 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.04,   0.04,   0.04,   0.04  }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.007,  0.007,  0.007,  0.007 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint80[6][8] = {
    {0.04,  0.04,  0.04,  0.04,  0.10,   0.10,   0.10,   0.10 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.03,   0.03,   0.03,   0.03 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.007,  0.007,  0.007,  0.007}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint70[6][8] = {
    {0.025, 0.025, 0.025, 0.025, 0.012,  0.012,  0.012,  0.012}, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.03,  0.03,  0.03,  0.03,  0.02,   0.02,   0.02,   0.02 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.005,  0.005,  0.005,  0.005}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint80NoHOverEE[6][8] = {
    {0.04,  0.04,  0.04,  0.04,  0.10,   0.10,   0.10,   0.10 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.03,   0.03,   0.03,   0.03 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.007,  0.007,  0.007,  0.007}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint70NoHOverEE[6][8] = {
    {0.025, 0.025, 0.025, 0.025, 0.10,   0.10,   0.10,   0.10 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.03,  0.03,  0.03,  0.03,  0.02,   0.02,   0.02,   0.02 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.005,  0.005,  0.005,  0.005}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  switch (idType) {
    case kCustomIdTight:    
      memcpy(fCuts,tightcuts,sizeof(fCuts));
      break;
    case kCustomIdLoose:
      memcpy(fCuts,loosecuts,sizeof(fCuts));
      break;
    case kVBTFWorkingPointFakeableId:
      memcpy(fCuts,VBTFWorkingPointFakeable,sizeof(fCuts));
      break;
    case kVetoId:
      memcpy(fCuts,Veto,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint95Id:
      memcpy(fCuts,VBTFWorkingPoint95,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint90Id:
      memcpy(fCuts,VBTFWorkingPoint90,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint85Id:
      memcpy(fCuts,VBTFWorkingPoint85,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint80Id:
      memcpy(fCuts,VBTFWorkingPoint80,sizeof(fCuts));
      break;
    case kVBTFWorkingPointLowPtId:
      if(ele->Pt() < 20)
        memcpy(fCuts,VBTFWorkingPoint70NoHOverEE,sizeof(fCuts));
      else
        memcpy(fCuts,VBTFWorkingPoint80NoHOverEE,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint70Id:
      memcpy(fCuts,VBTFWorkingPoint70,sizeof(fCuts));
      break;
    default:
      memset(fCuts,0,sizeof(fCuts));
      break;
  }

  // Based on RecoEgamma/ElectronIdentification/src/CutBasedElectronID.cc.
  Double_t eOverP = ele->ESuperClusterOverP();
  Double_t fBrem  = ele->FBrem();

  if ( (fCuts[5][0]>0.0) && (eOverP < fCuts[5][0]) && (fCuts[5][1]>0.0) && (fBrem < fCuts[5][1]))
    return kFALSE;

  if ( (fCuts[5][2]>0.0) && (eOverP < fCuts[5][2]*(1-fBrem)))
    return kFALSE;

  Int_t cat = 2;
  if ((ele->IsEB() && fBrem<0.06) || (ele->IsEE() && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  
  if(idType == kCustomIdLoose){
    double eleOneOverEMinusOneOverP = TMath::Abs((1.0/(ele->EcalEnergy())) - (1.0 / ele->P()));
    if(eleOneOverEMinusOneOverP >= 0.05) return kFALSE;
  }

  if(ele->SCluster() == 0)
    return kFALSE;
  Double_t eSeedOverPin = ele->ESeedClusterOverPIn();
  Double_t hOverE       = ele->HadronicOverEm();
  Double_t sigmaee      = ele->CoviEtaiEta();
  Double_t deltaPhiIn   = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn   = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());

  Int_t eb = 1;
  if (ele->IsEB()) 
    eb = 0;

  if (hOverE>fCuts[0][cat+4*eb])
    return kFALSE;

  if (sigmaee>fCuts[1][cat+4*eb])
    return kFALSE;

  if (deltaPhiIn>fCuts[2][cat+4*eb])
    return kFALSE; 

  if(deltaEtaIn>fCuts[3][cat+4*eb])
    return kFALSE; 

  if(eSeedOverPin<fCuts[4][cat+4*eb])
    return kFALSE;

  // Apply detector isolation at high pt only
  Bool_t isoCut = kTRUE;
  if(idType == kVBTFWorkingPointFakeableId){
    double isoEcal = ele->EcalRecHitIsoDr03();
    if(ele->IsEB()) isoEcal = isoEcal - 1.0;
    isoCut = (ele->TrackIsolationDr03() < ele->Pt()*0.2) &&
   	     (isoEcal                   < ele->Pt()*0.2) &&
   	     (ele->HcalTowerSumEtDr03() < ele->Pt()*0.2);
  }
  if(isoCut == kFALSE) return kFALSE;

  // Cuts only for pt<20 region and kVBTFWorkingPointLowPtId
  if(ele->Pt() < 20 && idType == kVBTFWorkingPointLowPtId) {
    Bool_t isGoodLowPtEl = fBrem > 0.15 ||
                          (ele->SCluster()->AbsEta() < 1.0 && eOverP > 0.95);
    if(!isGoodLowPtEl) return kFALSE;
  }
  
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassCustomIso(const Electron *ele, EElIsoType isoType,
                                    Bool_t useCombineIso) 
{
  Bool_t pass = kTRUE;
  Double_t fIsoCuts[4][2];          //!custom isolation cuts
  Double_t VBTFWorkingPoint95[4][2] = {
    {0.15 , 0.08   },   //TrkIso
    {2.00 , 0.06   },   //ECALIso
    {0.12 , 0.05   },   //HCALIso
    {0.15,  0.10   }   //Combined    
  };            

  Double_t VBTFWorkingPoint90[4][2] = {
    {0.12 , 0.05   },   //TrkIso
    {0.09 , 0.06   },   //ECALIso
    {0.10 , 0.03   },   //HCALIso
    {0.10,  0.07   }   //Combined    
  };            

  Double_t VBTFWorkingPoint85[4][2] = {
    {0.09 , 0.05   },   //TrkIso
    {0.08 , 0.05   },   //ECALIso
    {0.10 , 0.025  },   //HCALIso
    {0.09,  0.06   }   //Combined    
  };            

  Double_t VBTFWorkingPoint80[4][2] = {
    {0.09 , 0.04   },   //TrkIso
    {0.07 , 0.05   },   //ECALIso
    {0.10 , 0.025  },   //HCALIso
    {0.07,  0.06   }   //Combined    
  };            

  Double_t VBTFWorkingPoint70[4][2] = {
    {0.05 , 0.025  },   //TrkIso
    {0.06 , 0.025  },   //ECALIso
    {0.03 , 0.020  },   //HCALIso
    {0.04,  0.030  }   //Combined
  };            

  switch (isoType) {
    case kVBTFWorkingPoint95Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint95,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint90Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint90,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint85Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint85,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint80Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint80,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint70Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint70,sizeof(fIsoCuts));
      break;
    default:
      memset(fIsoCuts,0,sizeof(fIsoCuts));
      break;
  }

  Double_t trkIso  = ele->TrackIsolationDr03() / ele->Pt();
  Double_t ecalIso = ele->EcalRecHitIsoDr03() / ele->Pt();
  Double_t hcalIso = ele->HcalTowerSumEtDr03() / ele->Pt();
  Double_t combinedIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03();
  if(ele->IsEB()) combinedIso = ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03();
  combinedIso = combinedIso / ele->Pt();

  Int_t eb = 1;
  if (ele->IsEB()) 
    eb = 0;
 
  if(useCombineIso == kFALSE){
    if (trkIso>fIsoCuts[0][eb])
      pass = kFALSE;
    if (ecalIso>fIsoCuts[1][eb])
      pass = kFALSE;
    if (hcalIso>fIsoCuts[2][eb])
      pass = kFALSE;
  }
  else {
    if (combinedIso>fIsoCuts[3][eb])
      pass = kFALSE;
  }

  return pass;
}


//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassConversionFilter(const Electron *ele, 
                                           const DecayParticleCol *conversions,
                                           const BaseVertex *vtx,
                                           UInt_t nWrongHitsMax,
                                           Double_t probMin,
                                           Double_t lxyMin,
                                           Bool_t matchCkf,
                                           Bool_t requireArbitratedMerged,
                                           Double_t trkptMin) 
{
  Bool_t isGoodConversion = kFALSE;

  for (UInt_t ifc=0; ifc<conversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (conversions->At(ifc)->Daughter(d))->Trk();
      if (ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk) ) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (conversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || conversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (conversions->At(ifc)->LxyCorrected(vtx) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (conversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            if (trk->Pt()<trkptMin) isGoodConversion = kFALSE;
            const StableData *sd = dynamic_cast<const StableData*>
              (conversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }

    if (isGoodConversion == kTRUE) break;
    
  } // loop over all conversions 

  return !isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassConversionFilterPFAOD(const Electron *ele, 
                                           const DecayParticleCol *conversions,
                                           const BaseVertex *vtx,
                                           UInt_t nWrongHitsMax,
                                           Double_t probMin,
                                           Double_t lxyMin,
                                           Bool_t matchCkf,
                                           Bool_t requireArbitratedMerged,
                                           Double_t trkptMin) 
{

  Bool_t isGoodConversion = kFALSE;

  for (UInt_t ifc=0; ifc<conversions->GetEntries(); ifc++) {
    if(!(conversions->At(ifc)->Prob() > probMin) && (!requireArbitratedMerged || conversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) && (conversions->At(ifc)->LxyCorrected((BaseVertex*)vtx) > lxyMin)) continue;
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
      const ChargedParticle *pParticle = 0;
      pParticle = dynamic_cast<const ChargedParticle*>(conversions->At(ifc)->Daughter(d));
      if(pParticle == 0) continue;
      const Track* trk = 0;
      trk = pParticle->Trk();
      if(trk == 0) continue;
      if (ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk) ) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (conversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || conversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (conversions->At(ifc)->LxyCorrected(vtx) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
	  const ChargedParticle *pParticle = 0;
	  pParticle = dynamic_cast<const ChargedParticle*>(conversions->At(ifc)->Daughter(d));
	  if(pParticle == 0) continue;
	  const Track* trk = 0;
	  trk = pParticle->Trk();
	  if(trk == 0) continue;
          if (trk) {
            if (trk->Pt()<trkptMin) isGoodConversion = kFALSE;
            const StableData *sd = dynamic_cast<const StableData*>
              (conversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }
    
    if (isGoodConversion == kTRUE) break;
    
  } // loop over all conversions 
  
  return !isGoodConversion;
}
//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const VertexCol *vertices, Double_t fD0Cut, Int_t nVertex) 
{
  Bool_t d0cut = kFALSE;
  
  Double_t d0_real = 1e30;
  Double_t dz_real = 1e30;
  Int_t closestVtx = 0;

  if( nVertex >= (int) vertices->GetEntries() )
    nVertex = vertices->GetEntries() - 1;
  
  if(nVertex >= 0) d0_real = TMath::Abs(ele->GsfTrk()->D0Corrected(*vertices->At(nVertex)));
  else            {
    Double_t distVtx = 999.0;
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      double dz = TMath::Abs(ele->GsfTrk()->DzCorrected(*vertices->At(nv)));
      if(dz < distVtx) {
	distVtx    = dz;
        closestVtx = nv;
      }
      d0_real = TMath::Abs(ele->GsfTrk()->D0Corrected(*vertices->At(nv)));
      //printf("ming sync check nv:%d dz:%f d0:%f\n",nv,dz,d0_real);
    }
    d0_real = TMath::Abs(ele->GsfTrk()->D0Corrected(*vertices->At(closestVtx)));
    dz_real = TMath::Abs(ele->GsfTrk()->DzCorrected(*vertices->At(closestVtx)));
  }
  
  //printf("ming sync check electonidmod vtxind:%d vtx_x:%f vtx_y:%f vtx_z:%f dz:%f d0:%f\n",closestVtx,vertices->At(closestVtx)->X(),vertices->At(closestVtx)->Y(),vertices->At(closestVtx)->Z(),dz_real,d0_real);
  
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const BeamSpotCol *beamspots, Double_t fD0Cut) 
{
  Bool_t d0cut = kFALSE;
  // d0 cut
  Double_t d0_real = 99999;
  for(UInt_t i0 = 0; i0 < beamspots->GetEntries(); i0++) {
    Double_t pD0 = ele->GsfTrk()->D0Corrected(*beamspots->At(i0));
    if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassDZCut(const Electron *ele, const VertexCol *vertices, Double_t fDZCut, Int_t nVertex) 
{
  Bool_t dzcut = kFALSE;

  Double_t distVtx = 999.0;

  if( nVertex >= (int) vertices->GetEntries() )
    nVertex = vertices->GetEntries()-1;
  
  if(nVertex >= 0) distVtx = TMath::Abs(ele->GsfTrk()->DzCorrected(*vertices->At(nVertex)));
  else {
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      double dz = TMath::Abs(ele->GsfTrk()->DzCorrected(*vertices->At(nv)));
      if(dz < distVtx) {
        distVtx	 = dz;
      }
    }
  }
  
  if(distVtx < fDZCut) dzcut = kTRUE;
  
  return dzcut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassChargeFilter(const Electron *ele) 
{
  Bool_t passChargeFilter = kTRUE;
  if(ele->TrackerTrk() &&
     ele->TrackerTrk()->Charge() != ele->Charge()) passChargeFilter = kFALSE;

  return passChargeFilter;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassSpikeRemovalFilter(const Electron *ele) 
{
  Bool_t passSpikeRemovalFilter = kTRUE;
  if(ele->SCluster() &&
     ele->SCluster()->Seed()->Energy() > 5.0 && 
     ele->SCluster()->Seed()->EMax() / ele->SCluster()->Seed()->E3x3() > 0.95
    ) {
    passSpikeRemovalFilter = kFALSE;
  }

  // For Now Only use the EMax/E3x3 prescription.
  //   if(ele->SCluster()->Seed()->Energy() > 5.0 && 
  //      (1 - (ele->SCluster()->Seed()->E1x3() + ele->SCluster()->Seed()->E3x1() - 2*ele->SCluster()->Seed()->EMax())) > 0.95
  //     ) {
  //     passSpikeRemovalFilter = kFALSE;
  //   }
    
  return passSpikeRemovalFilter;
}

Bool_t ElectronTools::PassTriggerMatching(const Electron *ele, const TriggerObjectCol *trigobjs)
{
  
  for (UInt_t i=0; i<trigobjs->GetEntries(); ++i) {
    const TriggerObject *trigobj = trigobjs->At(i);
    if (trigobj->TriggerType()==TriggerObject::TriggerCluster || trigobj->TriggerType()==TriggerObject::TriggerElectron || trigobj->TriggerType()==TriggerObject::TriggerPhoton) {
      if (MathUtils::DeltaR(ele->SCluster(),trigobj)<0.3) {
        return kTRUE;
      }
    }
  }
  
  return kFALSE;
  
  
}

//--------------------------------------------------------------------------------------------------
Int_t ElectronTools::Classify(const Electron *ele) {
  
  double eta    = ele->AbsEta();
  double eOverP = ele->ESuperClusterOverP();
  double fBrem  = ele->FBrem();

  int cat = -1;
  if (ele->IsEB()) {
    if ((fBrem >= 0.12) and (eOverP > 0.9) and (eOverP < 1.2))
      cat = 0;
    else if (((eta >  .445   and eta <  .45  ) or
  	      (eta >  .79    and eta <  .81  ) or
  	      (eta > 1.137   and eta < 1.157 ) or
  	      (eta > 1.47285 and eta < 1.4744)))
      cat = 6;
    else if (ele->IsTrackerDriven() and !ele->IsEcalDriven())
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) and (eOverP > 0.82) and (eOverP < 1.22))
      cat = 3;
    else if (eta > 1.5 and eta <  1.58)
      cat = 7;
    else if (ele->IsTrackerDriven() and !ele->IsEcalDriven())
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else
      cat = 5;
  }

  return cat;
}

//--------------------------------------------------------------------------------------------------
Int_t ElectronTools::PassTightId(const Electron *ele, const VertexCol *vertices, 
                                 const DecayParticleCol *conversions, const Int_t typeCuts,
				 Double_t beta){

// original code on
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoEgamma/ElectronIdentification/src/CutBasedElectronID.cc
// beta must be computed with a DR cone of 0.4

  Double_t scEt   = ele->SCluster()->Et();
  Double_t scEta  = ele->SCluster()->Eta();

  Double_t fBrem = ele->FBrem();
  Double_t hOverE = ele->HadronicOverEm();
  Double_t sigmaee = ele->CoviEtaiEta();
  Double_t deltaPhiIn = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());
  Double_t eSeedOverPin = ele->ESeedClusterOverPIn(); 

  Int_t mishits = ele->BestTrk()->NExpectedHitsInner();
  Double_t tkIso   = ele->TrackIsolationDr03();
  Double_t ecalIso = ele->EcalRecHitIsoDr04()*beta;
  Double_t hcalIso  = ele->HcalTowerSumEtDr04()*beta;

  int cat = Classify(ele);
  int eb;

  if (ele->IsEB()) 
    eb = 0;
  else 
    eb = 1; 

  // Medium cuts
  Double_t cutdcotdistMedium[9] = {
  3.32e-02, 2.92e-02, 2.49e-02, 3.92e-02, 3.41e-02, 3.96e-02, 2.91e-02, 3.95e-02, 7.71e-03};
  Double_t cutdetainMedium[9] = {
  1.33e-02, 4.48e-03, 9.22e-03, 1.54e-02, 7.26e-03, 1.24e-02, 1.29e-02, 3.84e-02, 1.88e-02};
  Double_t cutdetainlMedium[9] = {
  1.21e-02, 4.22e-03, 9.18e-03, 1.61e-02, 6.45e-03, 1.16e-02, 1.23e-02, 6.20e-02, 2.43e-02};
  Double_t cutdphiinMedium[9] = {
  7.09e-02, 2.43e-01, 2.96e-01, 7.98e-02, 2.35e-01, 2.76e-01, 3.42e-01, 4.04e-01, 2.99e-01};
  Double_t cutdphiinlMedium[9] = {
  7.42e-02, 2.43e-01, 2.97e-01, 9.12e-02, 2.26e-01, 2.76e-01, 3.34e-01, 5.58e-01, 2.91e-01};
  Double_t cuteseedopcorMedium[9] = {
  6.42e-01, 9.44e-01, 4.53e-01, 7.62e-01, 3.67e-01, 5.57e-01, 1.98e-01, 9.15e-01, 6.28e-02};
  Double_t cutfmishitsMedium[9] = {
  4.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01};
  Double_t cuthoeMedium[9] = {
  1.96e-01, 6.30e-02, 1.48e-01, 3.66e-01, 5.66e-02, 1.45e-01, 4.29e-01, 4.28e-01, 3.99e-01};
  Double_t cuthoelMedium[9] = {
  2.19e-01, 6.19e-02, 1.47e-01, 3.58e-01, 4.61e-02, 1.46e-01, 3.26e-01, 3.81e-01, 3.89e-01};
  Double_t cutip_gsfMedium[9] = {
  2.45e-02, 9.74e-02, 1.48e-01, 5.49e-02, 5.65e-01, 3.33e-01, 2.04e-01, 5.41e-01, 1.21e-01};
  Double_t cutip_gsflMedium[9] = {
  1.92e-02, 9.81e-02, 1.33e-01, 4.34e-02, 5.65e-01, 3.24e-01, 2.33e-01, 4.30e-01, 6.44e-02};
  Double_t cutiso_sumMedium[9] = {
  1.44e+01, 1.12e+01, 1.09e+01, 1.08e+01, 6.35e+00, 9.78e+00, 1.30e+01, 1.62e+01, 1.96e+00};
  Double_t cutiso_sumoetMedium[9] = {
  1.01e+01, 6.41e+00, 6.00e+00, 8.14e+00, 3.90e+00, 4.76e+00, 6.86e+00, 6.48e+00, 1.74e+01};
  Double_t cutiso_sumoetlMedium[9] = {
  9.44e+00, 7.67e+00, 7.15e+00, 7.34e+00, 3.35e+00, 4.70e+00, 8.32e+00, 7.55e+00, 6.25e+00};
  Double_t cutseeMedium[9] = {
  1.30e-02, 1.09e-02, 1.18e-02, 3.94e-02, 3.04e-02, 3.28e-02, 1.00e-02, 3.73e-02, 6.69e-02};
  Double_t cutseelMedium[9] = {
  1.42e-02, 1.11e-02, 1.29e-02, 4.32e-02, 2.96e-02, 3.82e-02, 1.01e-02, 4.45e-02, 1.19e-01};

  // Tight cuts
  Double_t cutdcotdistTight[9] = {
  2.68e-02, 2.36e-02, 2.21e-02, 3.72e-02, 3.17e-02, 3.61e-02, 2.55e-02, 3.75e-02, 2.16e-04};
  Double_t cutdetainTight[9] = {
  8.92e-03, 3.96e-03, 8.50e-03, 1.34e-02, 6.27e-03, 1.05e-02, 1.12e-02, 3.09e-02, 1.88e-02};
  Double_t cutdetainlTight[9] = {
  9.23e-03, 3.77e-03, 8.70e-03, 1.39e-02, 5.60e-03, 9.40e-03, 1.07e-02, 6.20e-02, 4.10e-03};
  Double_t cutdphiinTight[9] = {
  6.37e-02, 1.53e-01, 2.90e-01, 7.69e-02, 1.81e-01, 2.34e-01, 3.42e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlTight[9] = {
  6.92e-02, 2.33e-01, 2.96e-01, 8.65e-02, 1.85e-01, 2.76e-01, 3.34e-01, 3.53e-01, 2.90e-01};
  Double_t cuteseedopcorTight[9] = {
  6.52e-01, 9.69e-01, 9.12e-01, 7.79e-01, 3.67e-01, 6.99e-01, 3.28e-01, 9.67e-01, 5.89e-01};
  Double_t cutfmishitsTight[9] = {
  4.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeTight[9] = {
  1.74e-01, 4.88e-02, 1.46e-01, 3.64e-01, 4.93e-02, 1.45e-01, 4.29e-01, 4.20e-01, 3.99e-01};
  Double_t cuthoelTight[9] = {
  2.19e-01, 5.25e-02, 1.47e-01, 3.57e-01, 4.25e-02, 1.45e-01, 3.26e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfTight[9] = {
  1.58e-02, 8.25e-02, 1.15e-01, 4.05e-02, 5.40e-01, 1.51e-01, 7.74e-02, 4.17e-01, 7.80e-02};
  Double_t cutip_gsflTight[9] = {
  1.27e-02, 6.26e-02, 9.68e-02, 3.02e-02, 5.65e-01, 1.46e-01, 7.90e-02, 4.10e-01, 4.79e-02};
  Double_t cutiso_sumTight[9] = {
  1.23e+01, 9.77e+00, 1.01e+01, 9.77e+00, 6.13e+00, 7.55e+00, 1.30e+01, 1.62e+01, 1.78e+00};
  Double_t cutiso_sumoetTight[9] = {
  7.75e+00, 5.45e+00, 5.67e+00, 5.97e+00, 3.17e+00, 3.86e+00, 6.06e+00, 5.31e+00, 1.05e+01};
  Double_t cutiso_sumoetlTight[9] = {
  7.56e+00, 5.08e+00, 5.77e+00, 5.74e+00, 2.37e+00, 3.32e+00, 4.97e+00, 5.46e+00, 3.82e+00};
  Double_t cutseeTight[9] = {
  1.16e-02, 1.07e-02, 1.08e-02, 3.49e-02, 2.89e-02, 3.08e-02, 9.87e-03, 3.37e-02, 4.40e-02};
  Double_t cutseelTight[9] = {
  1.27e-02, 1.08e-02, 1.13e-02, 4.19e-02, 2.81e-02, 3.02e-02, 9.76e-03, 4.28e-02, 2.98e-02};
  
  // SuperTight cuts
  Double_t cutdcotdistSuperTight[9] = {
  2.11e-02, 1.86e-02, 1.55e-02, 3.40e-02, 2.85e-02, 3.32e-02, 1.64e-02, 3.75e-02, 1.30e-04};
  Double_t cutdetainSuperTight[9] = {
  7.84e-03, 3.67e-03, 7.00e-03, 1.28e-02, 5.65e-03, 9.53e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlSuperTight[9] = {
  7.61e-03, 3.28e-03, 6.57e-03, 1.03e-02, 5.05e-03, 8.55e-03, 1.07e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinSuperTight[9] = {
  4.83e-02, 7.39e-02, 2.38e-01, 5.74e-02, 1.29e-01, 2.13e-01, 3.31e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlSuperTight[9] = {
  5.79e-02, 7.21e-02, 2.18e-01, 7.70e-02, 1.41e-01, 2.11e-01, 2.43e-01, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorSuperTight[9] = {
  7.32e-01, 9.77e-01, 9.83e-01, 8.55e-01, 4.31e-01, 7.35e-01, 4.18e-01, 9.99e-01, 5.89e-01};
  Double_t cutfmishitsSuperTight[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeSuperTight[9] = {
  9.19e-02, 4.11e-02, 1.42e-01, 3.35e-01, 3.82e-02, 1.41e-01, 4.29e-01, 4.01e-01, 3.99e-01};
  Double_t cuthoelSuperTight[9] = {
  7.51e-02, 3.81e-02, 1.41e-01, 3.32e-01, 3.10e-02, 1.43e-01, 2.35e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfSuperTight[9] = {
  1.42e-02, 2.66e-02, 1.06e-01, 3.38e-02, 3.23e-01, 1.07e-01, 7.74e-02, 2.32e-01, 7.80e-02};
  Double_t cutip_gsflSuperTight[9] = {
  1.15e-02, 2.72e-02, 8.41e-02, 2.49e-02, 4.17e-01, 1.02e-01, 7.90e-02, 1.69e-01, 4.79e-02};
  Double_t cutiso_sumSuperTight[9] = {
  8.95e+00, 8.18e+00, 8.75e+00, 7.47e+00, 5.43e+00, 5.87e+00, 8.16e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetSuperTight[9] = {
  6.45e+00, 5.14e+00, 4.99e+00, 5.21e+00, 2.65e+00, 3.12e+00, 4.52e+00, 4.72e+00, 3.68e+00};
  Double_t cutiso_sumoetlSuperTight[9] = {
  6.02e+00, 3.96e+00, 4.23e+00, 4.73e+00, 1.99e+00, 2.64e+00, 3.72e+00, 3.81e+00, 1.44e+00};
  Double_t cutseeSuperTight[9] = {
  1.09e-02, 1.05e-02, 1.05e-02, 3.24e-02, 2.81e-02, 2.95e-02, 9.77e-03, 2.75e-02, 2.95e-02};
  Double_t cutseelSuperTight[9] = {
  1.12e-02, 1.05e-02, 1.07e-02, 3.51e-02, 2.75e-02, 2.87e-02, 9.59e-03, 2.67e-02, 2.98e-02};

  // HyperTight1 cuts
  Double_t cutdcotdistHyperTight1[9] = {
  1.48e-02, 1.50e-02, 8.25e-03, 3.16e-02, 2.85e-02, 3.15e-02, 6.62e-03, 3.48e-02, 3.63e-06};
  Double_t cutdetainHyperTight1[9] = {
  6.51e-03, 3.51e-03, 5.53e-03, 9.16e-03, 5.30e-03, 8.28e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlHyperTight1[9] = {
  6.05e-03, 3.23e-03, 4.93e-03, 8.01e-03, 4.93e-03, 7.91e-03, 1.03e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinHyperTight1[9] = {
  4.83e-02, 4.91e-02, 2.30e-01, 3.48e-02, 7.44e-02, 2.04e-01, 9.95e-02, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlHyperTight1[9] = {
  4.74e-02, 4.51e-02, 2.18e-01, 2.99e-02, 7.37e-02, 2.11e-01, 9.99e-02, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorHyperTight1[9] = {
  7.72e-01, 9.90e-01, 1.01e+00, 8.55e-01, 9.11e-01, 7.72e-01, 9.17e-01, 1.06e+00, 7.63e-01};
  Double_t cutfmishitsHyperTight1[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeHyperTight1[9] = {
  6.17e-02, 3.70e-02, 1.41e-01, 2.91e-01, 3.82e-02, 1.34e-01, 4.19e-01, 3.87e-01, 3.93e-01};
  Double_t cuthoelHyperTight1[9] = {
  4.43e-02, 3.57e-02, 1.41e-01, 2.81e-01, 3.07e-02, 1.28e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight1[9] = {
  1.21e-02, 1.76e-02, 6.01e-02, 2.96e-02, 1.74e-01, 9.70e-02, 7.74e-02, 1.33e-01, 7.80e-02};
  Double_t cutip_gsflHyperTight1[9] = {
  1.01e-02, 1.56e-02, 6.87e-02, 2.13e-02, 1.25e-01, 8.16e-02, 7.90e-02, 1.30e-01, 4.79e-02};
  Double_t cutiso_sumHyperTight1[9] = {
  7.92e+00, 6.85e+00, 7.87e+00, 6.77e+00, 4.47e+00, 5.28e+00, 6.57e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight1[9] = {
  5.20e+00, 3.93e+00, 3.88e+00, 4.10e+00, 2.40e+00, 2.43e+00, 3.49e+00, 3.94e+00, 3.01e+00};
  Double_t cutiso_sumoetlHyperTight1[9] = {
  4.18e+00, 3.12e+00, 3.44e+00, 3.25e+00, 1.77e+00, 2.06e+00, 2.83e+00, 3.12e+00, 1.43e+00};
  Double_t cutseeHyperTight1[9] = {
  1.05e-02, 1.04e-02, 1.01e-02, 3.24e-02, 2.80e-02, 2.85e-02, 9.67e-03, 2.61e-02, 2.95e-02};
  Double_t cutseelHyperTight1[9] = {
  1.04e-02, 1.03e-02, 1.01e-02, 3.04e-02, 2.74e-02, 2.78e-02, 9.58e-03, 2.54e-02, 2.83e-02};

  // HyperTight2 cuts
  Double_t cutdcotdistHyperTight2[9] = {
  1.15e-02, 1.07e-02, 4.01e-03, 2.97e-02, 2.85e-02, 3.10e-02, 9.34e-04, 3.40e-02, 2.82e-07};
  Double_t cutdetainHyperTight2[9] = {
  5.29e-03, 2.56e-03, 4.89e-03, 7.89e-03, 5.30e-03, 7.37e-03, 8.91e-03, 9.36e-03, 5.94e-03};
  Double_t cutdetainlHyperTight2[9] = {
  4.48e-03, 2.59e-03, 4.42e-03, 6.54e-03, 4.93e-03, 6.98e-03, 8.49e-03, 9.06e-03, -4.81e-03};
  Double_t cutdphiinHyperTight2[9] = {
  2.41e-02, 3.83e-02, 1.48e-01, 2.91e-02, 3.15e-02, 1.57e-01, 8.90e-02, 1.02e-01, 2.81e-01};
  Double_t cutdphiinlHyperTight2[9] = {
  2.13e-02, 3.79e-02, 1.25e-01, 2.24e-02, 3.69e-02, 1.64e-01, 9.99e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight2[9] = {
  1.03e+00, 9.95e-01, 1.03e+00, 1.01e+00, 9.46e-01, 9.03e-01, 9.97e-01, 1.14e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight2[9] = {
  1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight2[9] = {
  4.94e-02, 3.45e-02, 1.40e-01, 2.02e-01, 3.82e-02, 1.19e-01, 1.23e-01, 3.82e-01, 2.50e-01};
  Double_t cuthoelHyperTight2[9] = {
  4.04e-02, 3.42e-02, 1.31e-01, 1.85e-01, 3.01e-02, 1.27e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight2[9] = {
  1.14e-02, 1.38e-02, 5.29e-02, 1.87e-02, 1.31e-01, 8.63e-02, 7.74e-02, 1.04e-01, 2.42e-02};
  Double_t cutip_gsflHyperTight2[9] = {
  9.83e-03, 1.35e-02, 4.27e-02, 1.72e-02, 1.25e-01, 7.92e-02, 7.90e-02, 1.30e-01, 3.40e-02};
  Double_t cutiso_sumHyperTight2[9] = {
  6.40e+00, 5.77e+00, 6.54e+00, 5.22e+00, 3.86e+00, 4.63e+00, 6.31e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight2[9] = {
  4.03e+00, 3.03e+00, 3.24e+00, 3.13e+00, 2.05e+00, 2.01e+00, 2.99e+00, 3.44e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight2[9] = {
  3.08e+00, 2.31e+00, 2.84e+00, 2.53e+00, 1.65e+00, 1.72e+00, 2.34e+00, 3.11e+00, 1.35e+00};
  Double_t cutseeHyperTight2[9] = {
  1.03e-02, 1.03e-02, 9.88e-03, 3.03e-02, 2.79e-02, 2.79e-02, 9.67e-03, 2.52e-02, 2.58e-02};
  Double_t cutseelHyperTight2[9] = {
  1.02e-02, 1.02e-02, 9.80e-03, 2.90e-02, 2.74e-02, 2.75e-02, 9.58e-03, 2.49e-02, 2.50e-02};

  // HyperTight3 cuts
  Double_t cutdcotdistHyperTight3[9] = {
  9.63e-03, 5.11e-03, 1.95e-04, 2.97e-02, 2.85e-02, 2.18e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight3[9] = {
  4.86e-03, 2.29e-03, 4.40e-03, 7.79e-03, 4.07e-03, 6.33e-03, 7.70e-03, 7.93e-03, 5.94e-03};
  Double_t cutdetainlHyperTight3[9] = {
  4.48e-03, 2.30e-03, 4.14e-03, 6.04e-03, 3.87e-03, 6.09e-03, 7.97e-03, 8.04e-03, -4.81e-03};
  Double_t cutdphiinHyperTight3[9] = {
  2.41e-02, 2.88e-02, 7.39e-02, 2.91e-02, 1.91e-02, 1.14e-01, 3.61e-02, 8.92e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight3[9] = {
  1.95e-02, 3.42e-02, 8.06e-02, 2.22e-02, 2.26e-02, 9.73e-02, 4.51e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight3[9] = {
  1.07e+00, 1.01e+00, 1.08e+00, 1.01e+00, 9.69e-01, 9.10e-01, 1.04e+00, 1.20e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight3[9] = {
  5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight3[9] = {
  3.52e-02, 3.45e-02, 1.33e-01, 1.88e-01, 2.72e-02, 1.19e-01, 9.28e-02, 2.46e-01, 2.50e-01};
  Double_t cuthoelHyperTight3[9] = {
  4.04e-02, 3.40e-02, 1.31e-01, 1.84e-01, 2.64e-02, 1.18e-01, 9.76e-02, 2.53e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight3[9] = {
  1.14e-02, 1.26e-02, 3.79e-02, 1.68e-02, 1.21e-01, 5.29e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight3[9] = {
  9.83e-03, 1.18e-02, 3.59e-02, 1.56e-02, 1.20e-01, 5.36e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight3[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.87e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight3[9] = {
  3.03e+00, 2.50e+00, 2.58e+00, 2.44e+00, 1.91e+00, 1.76e+00, 2.92e+00, 3.13e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight3[9] = {
  2.36e+00, 2.02e+00, 2.29e+00, 1.89e+00, 1.65e+00, 1.69e+00, 2.03e+00, 2.79e+00, 1.35e+00};
  Double_t cutseeHyperTight3[9] = {
  1.03e-02, 1.01e-02, 9.84e-03, 2.89e-02, 2.74e-02, 2.73e-02, 9.47e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight3[9] = {
  1.02e-02, 1.00e-02, 9.73e-03, 2.79e-02, 2.73e-02, 2.69e-02, 9.40e-03, 2.46e-02, 2.50e-02};

  // HyperTight4 cuts
  Double_t cutdcotdistHyperTight4[9] = {
  2.70e-04, 1.43e-04, 1.95e-04, 2.64e-03, 2.82e-02, 1.64e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight4[9] = {
  2.44e-03, 1.67e-03, 2.26e-03, 3.43e-03, 3.51e-03, 3.52e-03, 2.98e-03, 4.79e-03, 5.94e-03};
  Double_t cutdetainlHyperTight4[9] = {
  2.34e-03, 1.29e-03, 2.30e-03, 3.30e-03, 3.61e-03, 3.84e-03, 2.53e-03, 3.66e-03, -4.81e-03};
  Double_t cutdphiinHyperTight4[9] = {
  8.44e-03, 5.21e-03, 2.18e-02, 1.39e-02, 7.82e-03, 1.52e-02, 2.59e-02, 3.87e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight4[9] = {
  5.77e-03, 3.20e-03, 2.85e-02, 2.22e-02, 7.00e-03, 1.84e-02, 2.91e-02, 4.40e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight4[9] = {
  1.15e+00, 1.01e+00, 1.21e+00, 1.07e+00, 9.69e-01, 9.10e-01, 1.08e+00, 1.36e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight4[9] = {
  5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight4[9] = {
  2.39e-02, 2.68e-02, 2.12e-02, 1.03e-01, 9.92e-03, 7.07e-02, 7.12e-02, 1.48e-01, 2.50e-01};
  Double_t cuthoelHyperTight4[9] = {
  2.87e-02, 1.94e-02, 2.16e-02, 5.68e-02, 1.35e-02, 4.04e-02, 7.98e-02, 1.50e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight4[9] = {
  7.61e-03, 5.22e-03, 3.79e-02, 1.02e-02, 4.62e-02, 1.82e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight4[9] = {
  7.81e-03, 4.25e-03, 3.08e-02, 1.04e-02, 2.35e-02, 2.45e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight4[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.86e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight4[9] = {
  2.53e+00, 2.10e+00, 1.87e+00, 1.84e+00, 1.79e+00, 1.61e+00, 2.53e+00, 1.98e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight4[9] = {
  2.28e+00, 2.02e+00, 2.04e+00, 1.69e+00, 1.65e+00, 1.61e+00, 2.03e+00, 1.82e+00, 1.35e+00};
  Double_t cutseeHyperTight4[9] = {
  9.99e-03, 9.61e-03, 9.65e-03, 2.75e-02, 2.61e-02, 2.64e-02, 9.18e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight4[9] = {
  9.66e-03, 9.69e-03, 9.58e-03, 2.73e-02, 2.66e-02, 2.66e-02, 8.64e-03, 2.46e-02, 2.50e-02};

  Double_t cutdcotdist[9];
  Double_t cutdetain[9];
  Double_t cutdetainl[9];
  Double_t cutdphiin[9];
  Double_t cutdphiinl[9];
  Double_t cuteseedopcor[9];
  Double_t cutfmishits[9];
  Double_t cuthoe[9];
  Double_t cuthoel[9];
  Double_t cutip_gsf[9];
  Double_t cutip_gsfl[9];
  Double_t cutiso_sum[9];
  Double_t cutiso_sumoet[9];
  Double_t cutiso_sumoetl[9];
  Double_t cutsee[9];
  Double_t cutseel[9];
  if	 (typeCuts == 0) {
    memcpy(cutdcotdist   ,cutdcotdistMedium   ,sizeof(cutdcotdistMedium));
    memcpy(cutdetain     ,cutdetainMedium     ,sizeof(cutdetainMedium));
    memcpy(cutdetainl    ,cutdetainlMedium    ,sizeof(cutdetainlMedium));
    memcpy(cutdphiin     ,cutdphiinMedium     ,sizeof(cutdphiinMedium));
    memcpy(cutdphiinl    ,cutdphiinlMedium    ,sizeof(cutdphiinlMedium));
    memcpy(cuteseedopcor ,cuteseedopcorMedium ,sizeof(cuteseedopcorMedium));
    memcpy(cutfmishits   ,cutfmishitsMedium   ,sizeof(cutfmishitsMedium));
    memcpy(cuthoe        ,cuthoeMedium	      ,sizeof(cuthoeMedium));
    memcpy(cuthoel       ,cuthoelMedium	      ,sizeof(cuthoelMedium));
    memcpy(cutip_gsf     ,cutip_gsfMedium     ,sizeof(cutip_gsfMedium));
    memcpy(cutip_gsfl    ,cutip_gsflMedium    ,sizeof(cutip_gsflMedium));
    memcpy(cutiso_sum    ,cutiso_sumMedium    ,sizeof(cutiso_sumMedium));
    memcpy(cutiso_sumoet ,cutiso_sumoetMedium ,sizeof(cutiso_sumoetMedium));
    memcpy(cutiso_sumoetl,cutiso_sumoetlMedium,sizeof(cutiso_sumoetlMedium));
    memcpy(cutsee        ,cutseeMedium	      ,sizeof(cutseeMedium));
    memcpy(cutseel       ,cutseelMedium	      ,sizeof(cutseelMedium));
  }
  else if(typeCuts == 1) {
    memcpy(cutdcotdist   ,cutdcotdistTight   ,sizeof(cutdcotdistTight));
    memcpy(cutdetain     ,cutdetainTight     ,sizeof(cutdetainTight));
    memcpy(cutdetainl    ,cutdetainlTight    ,sizeof(cutdetainlTight));
    memcpy(cutdphiin     ,cutdphiinTight     ,sizeof(cutdphiinTight));
    memcpy(cutdphiinl    ,cutdphiinlTight    ,sizeof(cutdphiinlTight));
    memcpy(cuteseedopcor ,cuteseedopcorTight ,sizeof(cuteseedopcorTight));
    memcpy(cutfmishits   ,cutfmishitsTight   ,sizeof(cutfmishitsTight));
    memcpy(cuthoe        ,cuthoeTight	     ,sizeof(cuthoeTight));
    memcpy(cuthoel       ,cuthoelTight	     ,sizeof(cuthoelTight));
    memcpy(cutip_gsf     ,cutip_gsfTight     ,sizeof(cutip_gsfTight));
    memcpy(cutip_gsfl    ,cutip_gsflTight    ,sizeof(cutip_gsflTight));
    memcpy(cutiso_sum    ,cutiso_sumTight    ,sizeof(cutiso_sumTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetTight ,sizeof(cutiso_sumoetTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlTight,sizeof(cutiso_sumoetlTight));
    memcpy(cutsee        ,cutseeTight	     ,sizeof(cutseeTight));
    memcpy(cutseel       ,cutseelTight	     ,sizeof(cutseelTight));
  }
  else if(typeCuts == 2) {
    memcpy(cutdcotdist   ,cutdcotdistSuperTight   ,sizeof(cutdcotdistSuperTight));
    memcpy(cutdetain     ,cutdetainSuperTight     ,sizeof(cutdetainSuperTight));
    memcpy(cutdetainl    ,cutdetainlSuperTight    ,sizeof(cutdetainlSuperTight));
    memcpy(cutdphiin     ,cutdphiinSuperTight     ,sizeof(cutdphiinSuperTight));
    memcpy(cutdphiinl    ,cutdphiinlSuperTight    ,sizeof(cutdphiinlSuperTight));
    memcpy(cuteseedopcor ,cuteseedopcorSuperTight ,sizeof(cuteseedopcorSuperTight));
    memcpy(cutfmishits   ,cutfmishitsSuperTight   ,sizeof(cutfmishitsSuperTight));
    memcpy(cuthoe        ,cuthoeSuperTight	  ,sizeof(cuthoeSuperTight));
    memcpy(cuthoel       ,cuthoelSuperTight	  ,sizeof(cuthoelSuperTight));
    memcpy(cutip_gsf     ,cutip_gsfSuperTight     ,sizeof(cutip_gsfSuperTight));
    memcpy(cutip_gsfl    ,cutip_gsflSuperTight    ,sizeof(cutip_gsflSuperTight));
    memcpy(cutiso_sum    ,cutiso_sumSuperTight    ,sizeof(cutiso_sumSuperTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetSuperTight ,sizeof(cutiso_sumoetSuperTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlSuperTight,sizeof(cutiso_sumoetlSuperTight));
    memcpy(cutsee        ,cutseeSuperTight	  ,sizeof(cutseeSuperTight));
    memcpy(cutseel       ,cutseelSuperTight	  ,sizeof(cutseelSuperTight));
  }
  else if(typeCuts == 3) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight1   ,sizeof(cutdcotdistHyperTight1));
    memcpy(cutdetain     ,cutdetainHyperTight1     ,sizeof(cutdetainHyperTight1));
    memcpy(cutdetainl    ,cutdetainlHyperTight1    ,sizeof(cutdetainlHyperTight1));
    memcpy(cutdphiin     ,cutdphiinHyperTight1     ,sizeof(cutdphiinHyperTight1));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight1    ,sizeof(cutdphiinlHyperTight1));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight1 ,sizeof(cuteseedopcorHyperTight1));
    memcpy(cutfmishits   ,cutfmishitsHyperTight1   ,sizeof(cutfmishitsHyperTight1));
    memcpy(cuthoe        ,cuthoeHyperTight1	  ,sizeof(cuthoeHyperTight1));
    memcpy(cuthoel       ,cuthoelHyperTight1	  ,sizeof(cuthoelHyperTight1));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight1     ,sizeof(cutip_gsfHyperTight1));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight1    ,sizeof(cutip_gsflHyperTight1));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight1    ,sizeof(cutiso_sumHyperTight1));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight1 ,sizeof(cutiso_sumoetHyperTight1));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight1,sizeof(cutiso_sumoetlHyperTight1));
    memcpy(cutsee        ,cutseeHyperTight1	  ,sizeof(cutseeHyperTight1));
    memcpy(cutseel       ,cutseelHyperTight1	  ,sizeof(cutseelHyperTight1));
  }
  else if(typeCuts == 4) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight2   ,sizeof(cutdcotdistHyperTight2));
    memcpy(cutdetain     ,cutdetainHyperTight2     ,sizeof(cutdetainHyperTight2));
    memcpy(cutdetainl    ,cutdetainlHyperTight2    ,sizeof(cutdetainlHyperTight2));
    memcpy(cutdphiin     ,cutdphiinHyperTight2     ,sizeof(cutdphiinHyperTight2));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight2    ,sizeof(cutdphiinlHyperTight2));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight2 ,sizeof(cuteseedopcorHyperTight2));
    memcpy(cutfmishits   ,cutfmishitsHyperTight2   ,sizeof(cutfmishitsHyperTight2));
    memcpy(cuthoe        ,cuthoeHyperTight2	  ,sizeof(cuthoeHyperTight2));
    memcpy(cuthoel       ,cuthoelHyperTight2	  ,sizeof(cuthoelHyperTight2));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight2     ,sizeof(cutip_gsfHyperTight2));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight2    ,sizeof(cutip_gsflHyperTight2));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight2    ,sizeof(cutiso_sumHyperTight2));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight2 ,sizeof(cutiso_sumoetHyperTight2));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight2,sizeof(cutiso_sumoetlHyperTight2));
    memcpy(cutsee        ,cutseeHyperTight2	  ,sizeof(cutseeHyperTight2));
    memcpy(cutseel       ,cutseelHyperTight2	  ,sizeof(cutseelHyperTight2));
  }
  else if(typeCuts == 5) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight3   ,sizeof(cutdcotdistHyperTight3));
    memcpy(cutdetain     ,cutdetainHyperTight3     ,sizeof(cutdetainHyperTight3));
    memcpy(cutdetainl    ,cutdetainlHyperTight3    ,sizeof(cutdetainlHyperTight3));
    memcpy(cutdphiin     ,cutdphiinHyperTight3     ,sizeof(cutdphiinHyperTight3));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight3    ,sizeof(cutdphiinlHyperTight3));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight3 ,sizeof(cuteseedopcorHyperTight3));
    memcpy(cutfmishits   ,cutfmishitsHyperTight3   ,sizeof(cutfmishitsHyperTight3));
    memcpy(cuthoe        ,cuthoeHyperTight3	  ,sizeof(cuthoeHyperTight3));
    memcpy(cuthoel       ,cuthoelHyperTight3	  ,sizeof(cuthoelHyperTight3));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight3     ,sizeof(cutip_gsfHyperTight3));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight3    ,sizeof(cutip_gsflHyperTight3));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight3    ,sizeof(cutiso_sumHyperTight3));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight3 ,sizeof(cutiso_sumoetHyperTight3));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight3,sizeof(cutiso_sumoetlHyperTight3));
    memcpy(cutsee        ,cutseeHyperTight3	  ,sizeof(cutseeHyperTight3));
    memcpy(cutseel       ,cutseelHyperTight3	  ,sizeof(cutseelHyperTight3));
  }
  else if(typeCuts == 6) {
    memcpy(cutdcotdist   ,cutdcotdistHyperTight4   ,sizeof(cutdcotdistHyperTight4));
    memcpy(cutdetain     ,cutdetainHyperTight4     ,sizeof(cutdetainHyperTight4));
    memcpy(cutdetainl    ,cutdetainlHyperTight4    ,sizeof(cutdetainlHyperTight4));
    memcpy(cutdphiin     ,cutdphiinHyperTight4     ,sizeof(cutdphiinHyperTight4));
    memcpy(cutdphiinl    ,cutdphiinlHyperTight4    ,sizeof(cutdphiinlHyperTight4));
    memcpy(cuteseedopcor ,cuteseedopcorHyperTight4 ,sizeof(cuteseedopcorHyperTight4));
    memcpy(cutfmishits   ,cutfmishitsHyperTight4   ,sizeof(cutfmishitsHyperTight4));
    memcpy(cuthoe        ,cuthoeHyperTight4	  ,sizeof(cuthoeHyperTight4));
    memcpy(cuthoel       ,cuthoelHyperTight4	  ,sizeof(cuthoelHyperTight4));
    memcpy(cutip_gsf     ,cutip_gsfHyperTight4     ,sizeof(cutip_gsfHyperTight4));
    memcpy(cutip_gsfl    ,cutip_gsflHyperTight4    ,sizeof(cutip_gsflHyperTight4));
    memcpy(cutiso_sum    ,cutiso_sumHyperTight4    ,sizeof(cutiso_sumHyperTight4));
    memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight4 ,sizeof(cutiso_sumoetHyperTight4));
    memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight4,sizeof(cutiso_sumoetlHyperTight4));
    memcpy(cutsee        ,cutseeHyperTight4	  ,sizeof(cutseeHyperTight4));
    memcpy(cutseel       ,cutseelHyperTight4	  ,sizeof(cutseelHyperTight4));
  }
  else {
    return 0;
  }
  int result = 0;
  
  const int ncuts = 10;
  std::vector<bool> cut_results(ncuts, false);
  
  float iso_sum = tkIso + ecalIso + hcalIso;
  if(fabs(scEta)>1.5) 
    iso_sum += (fabs(scEta)-1.5)*1.09;
  
  float iso_sumoet = iso_sum*(40./scEt);
  
  float eseedopincor = eSeedOverPin + fBrem;
  if(fBrem < 0)
    eseedopincor = eSeedOverPin;

  float dist = (TMath::Abs(ele->ConvPartnerDist())      == -9999.? 9999:TMath::Abs(ele->ConvPartnerDist()));
  float dcot = (TMath::Abs(ele->ConvPartnerDCotTheta()) == -9999.? 9999:TMath::Abs(ele->ConvPartnerDCotTheta()));

  float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);

  Double_t ip = 99999;
  for(UInt_t i0 = 0; i0 < vertices->GetEntries(); i0++) {
    if(vertices->At(i0)->NTracks() > 0){
      Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
      ip = TMath::Abs(pD0);
      break;
    }
  }

  for (int cut=0; cut<ncuts; cut++) {
    switch (cut) {
    case 0:
      cut_results[cut] = compute_cut(fabs(deltaEtaIn), scEt, cutdetainl[cat], cutdetain[cat]);
      break;
    case 1:
      cut_results[cut] = compute_cut(fabs(deltaPhiIn), scEt, cutdphiinl[cat], cutdphiin[cat]);
      break;
    case 2:
      cut_results[cut] = (eseedopincor > cuteseedopcor[cat]);
      break;
    case 3:
      cut_results[cut] = compute_cut(hOverE, scEt, cuthoel[cat], cuthoe[cat]);
      break;
    case 4:
      cut_results[cut] = compute_cut(sigmaee, scEt, cutseel[cat], cutsee[cat]);
      break;
    case 5:
      cut_results[cut] = compute_cut(iso_sumoet, scEt, cutiso_sumoetl[cat], cutiso_sumoet[cat]);
      break;
    case 6:
      cut_results[cut] = (iso_sum < cutiso_sum[cat]);
      break;
    case 7:
      cut_results[cut] = compute_cut(fabs(ip), scEt, cutip_gsfl[cat], cutip_gsf[cat]);
      break;
    case 8:
      cut_results[cut] = (mishits < cutfmishits[cat]);
      break;
    case 9:
      cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
      break;
    }
  }
  
  // ID part
  if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4])
    result = result + 1;
  
  // ISO part
  if (cut_results[5] & cut_results[6])
    result = result + 2;
  
  // IP part
  if (cut_results[7])
    result = result + 8;
  
  // Conversion part
  if (cut_results[8] and cut_results[9])
    result = result + 4;
  
  return result;
}

//--------------------------------------------------------------------------------------------------
bool ElectronTools::compute_cut(double x, double et, double cut_min, double cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronTools::Likelihood(ElectronLikelihood *LH, const Electron *ele) 
{
  if (!LH) {
    std::cout << "Error: Likelihood not properly initialized\n"; 
    return -9999;
  }

  LikelihoodMeasurements measurements;
  measurements.pt = ele->Pt();
  if (ele->IsEB() && ele->AbsEta()<1.0) measurements.subdet = 0;
  else if (ele->IsEB())                 measurements.subdet = 1;
  else                                  measurements.subdet = 2;
  measurements.deltaPhi = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  measurements.deltaEta = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());
  measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout();
  measurements.eSuperClusterOverP = ele->ESuperClusterOverP();
  measurements.hadronicOverEm = ele->HadronicOverEm();
  measurements.sigmaIEtaIEta = ele->CoviEtaiEta();
  measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi());
  measurements.fBrem = ele->FBrem();
  measurements.nBremClusters = ele->NumberOfClusters() - 1;
  //measurements.OneOverEMinusOneOverP = (1.0 / ele->SCluster()->Energy()) - (1.0 / ele->BestTrk()->P());
  measurements.OneOverEMinusOneOverP = (1.0 / ele->ESuperClusterOverP() / ele->BestTrk()->P()) - (1.0 / ele->BestTrk()->P());
  double likelihood = LH->result(measurements);

  double newLik = 0.0;
  if     (likelihood<=0) newLik = -20.0;
  else if(likelihood>=1) newLik =  20.0;
  else                   newLik = log(likelihood/(1.0-likelihood));

  Bool_t isDebug = kFALSE;
  if(isDebug == kTRUE){
    printf("LIKELIHOOD: %f %d %f %f %f %f %f %f %f %f %d %f %f %f - %f %f\n",measurements.pt,measurements.subdet,
    measurements.deltaPhi          ,measurements.deltaEta      ,measurements.eSeedClusterOverPout,
    measurements.eSuperClusterOverP,measurements.hadronicOverEm,measurements.sigmaIEtaIEta,
    measurements.sigmaIPhiIPhi     ,measurements.fBrem         ,measurements.nBremClusters,
    measurements.OneOverEMinusOneOverP,ele->SCluster()->Energy(),ele->BestTrk()->P(),
    likelihood,newLik);
  }

  return newLik;

}

Double_t ElectronTools::ElectronEffectiveArea(EElectronEffectiveAreaType type, Double_t SCEta, EElectronEffectiveAreaTarget EffectiveAreaTarget) {

  Double_t EffectiveArea = 0;
  
  if(type == ElectronTools::kEleNeutralIso04){
    if     (fabs(SCEta) < 1.0                           ) EffectiveArea = 0.208;
    else if(fabs(SCEta) >= 1.0   && fabs(SCEta) < 1.479 ) EffectiveArea = 0.209;
    else if(fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.000 ) EffectiveArea = 0.115;
    else if(fabs(SCEta) >= 2.000 && fabs(SCEta) < 2.200 ) EffectiveArea = 0.143;
    else if(fabs(SCEta) >= 2.200 && fabs(SCEta) < 2.300 ) EffectiveArea = 0.183;
    else if(fabs(SCEta) >= 2.300 && fabs(SCEta) < 2.400 ) EffectiveArea = 0.194;
    else                                                  EffectiveArea = 0.261;
    return EffectiveArea;
  }

  if (fabs(SCEta) < 1.0) {
    if (type == ElectronTools::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso03) EffectiveArea = 0.017;
    if (type == ElectronTools::kEleGammaIso03) EffectiveArea = 0.045;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.014;
    if (type == ElectronTools::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso04) EffectiveArea = 0.034;
    if (type == ElectronTools::kEleGammaIso04) EffectiveArea = 0.079;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.014;
    if (type == ElectronTools::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleHoverE) EffectiveArea = 0.00016;
    if (type == ElectronTools::kEleHcalDepth1OverEcal) EffectiveArea = 0.00016;
    if (type == ElectronTools::kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) {
    if (type == ElectronTools::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso03) EffectiveArea = 0.025;
    if (type == ElectronTools::kEleGammaIso03) EffectiveArea = 0.052;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.030;
    if (type == ElectronTools::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso04) EffectiveArea = 0.050;
    if (type == ElectronTools::kEleGammaIso04) EffectiveArea = 0.073;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.030;
    if (type == ElectronTools::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleHoverE) EffectiveArea = 0.00022;
    if (type == ElectronTools::kEleHcalDepth1OverEcal) EffectiveArea = 0.00022;
    if (type == ElectronTools::kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) {
    if (type == ElectronTools::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso03) EffectiveArea = 0.030;
    if (type == ElectronTools::kEleGammaIso03) EffectiveArea = 0.170;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.134;
    if (type == ElectronTools::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso04) EffectiveArea = 0.060;
    if (type == ElectronTools::kEleGammaIso04) EffectiveArea = 0.187;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.134;
    if (type == ElectronTools::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleHoverE) EffectiveArea = 0.00030;
    if (type == ElectronTools::kEleHcalDepth1OverEcal) EffectiveArea = 0.00026;
    if (type == ElectronTools::kEleHcalDepth2OverEcal) EffectiveArea = 0.00002;        
  } else if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.25 ) {
    if (type == ElectronTools::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso03) EffectiveArea = 0.022;
    if (type == ElectronTools::kEleGammaIso03) EffectiveArea = 0.623;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.516;
    if (type == ElectronTools::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso04) EffectiveArea = 0.055;
    if (type == ElectronTools::kEleGammaIso04) EffectiveArea = 0.659;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.517;
    if (type == ElectronTools::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleHoverE) EffectiveArea = 0.00054;
    if (type == ElectronTools::kEleHcalDepth1OverEcal) EffectiveArea = 0.00045;
    if (type == ElectronTools::kEleHcalDepth2OverEcal) EffectiveArea = 0.00003;
  } else if (fabs(SCEta) >= 2.25 && fabs(SCEta) < 2.5 ) {
    if (type == ElectronTools::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso03) EffectiveArea = 0.018;
    if (type == ElectronTools::kEleGammaIso03) EffectiveArea = 1.198;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 1.049;
    if (type == ElectronTools::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleNeutralHadronIso04) EffectiveArea = 0.073;
    if (type == ElectronTools::kEleGammaIso04) EffectiveArea = 1.258;
    if (type == ElectronTools::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 1.051;
    if (type == ElectronTools::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronTools::kEleHoverE) EffectiveArea = 0.00082;
    if (type == ElectronTools::kEleHcalDepth1OverEcal) EffectiveArea = 0.00066;
    if (type == ElectronTools::kEleHcalDepth2OverEcal) EffectiveArea = 0.00004;
  }
    
  //NoCorrections
  if (EffectiveAreaTarget == kEleEANoCorr) {
    return 0.0;
  }
  else if (EffectiveAreaTarget == kEleEAData2012) {
    if (type == kEleGammaIso03) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.122;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.147;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.106;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.138;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.221;
      if (fabs(SCEta) >= 2.4 ) EffectiveArea = 0.211;
    }
    if (type == kEleGammaIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.176;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.206;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.094;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.172;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.244;
      if (fabs(SCEta) >= 2.4 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.333;
      if (fabs(SCEta) >= 2.4 ) EffectiveArea = 0.348;
    }
    if (type == kEleNeutralHadronIso03) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.4 ) EffectiveArea = 0.019;
    }
    if (type == kEleNeutralHadronIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.022;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.052;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.063;
      if (fabs(SCEta) >= 2.4 ) EffectiveArea = 0.028;
    }
    if (type == kEleGammaAndNeutralHadronIso03) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.130;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.137;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.067;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.089;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.107;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.110;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.138;
    }				
    if (type == kEleGammaAndNeutralHadronIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.208;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.209;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.115;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.143;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.183;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.194;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.261;
    }				
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.013;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.052;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.102;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.034;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.042;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.053;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.069;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.033;
     if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.056;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.065;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.074;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.008;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.010;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.028;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.063;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.143;
    }
  } 

  //2011 Data Effective Areas
  else if (EffectiveAreaTarget == kEleEAData2011) {
    if (type == kEleGammaAndNeutralHadronIso03) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.100;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.120;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.085;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.110;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.120;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.120;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.130;
    }
    if (type == kEleGammaAndNeutralHadronIso04) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.180;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.200;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.150;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.190;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.210;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.220;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.290;
    }
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.041;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.048;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.029;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.054;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.028;
     if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.059;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.003;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.019;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.025;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.084;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.107;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.228;
    }
  } 

  //Summer11 MC Effective Areas
  else if (EffectiveAreaTarget == kEleEASummer11MC) {
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.023;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.046;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.046;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.059;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.058;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.073;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.083;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.053;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.048;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.062;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.085;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.118;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.004;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.021;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.040;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.086;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.105;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.169;
    }
  } 
  
  //Fall11 MC Effective Areas
  else if (EffectiveAreaTarget == kEleEAFall11MC) {
    if (type == kEleGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.012;
    }
    if (type == kEleGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.068;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.067;
    }
    if (type == kEleGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.080;
    }
    if (type == kEleGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.040;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.083;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.123;
    }
    if (type == kEleGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.095;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.123;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.133;
    }
    if (type == kEleNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.000;
    }
    if (type == kEleNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.007;
    }
    if (type == kEleNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.021;
    }
    if (type == kEleNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.110;
    }
    if (type == kEleNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.123;
      if (fabs(SCEta) >= 2.4) EffectiveArea = 0.220;
    }
  }

  return EffectiveArea;  
}


Bool_t ElectronTools::PassHggLeptonTagID(const Electron* ele) {
  
  float dist = ( ele->ConvPartnerDist()      == -9999.? 9999:TMath::Abs(ele->ConvPartnerDist()));
  float dcot = ( ele->ConvPartnerDCotTheta() == -9999.? 9999:TMath::Abs(ele->ConvPartnerDCotTheta()));  
  
  if (dist < 0.02) return false;
  if (dcot < 0.02) return false;
  
  int numInnerHits = ele->Trk()->NExpectedHitsInner();
  if( numInnerHits > 1 ) return false;

  float coviEtaiEta = ele->CoviEtaiEta();
  if( ele->SCluster()->AbsEta() < 1.5 && coviEtaiEta > 0.01 ) return false;
  else if( ele->SCluster()->AbsEta() > 1.5 && coviEtaiEta > 0.031 ) return false;

  Double_t deltaPhiIn   = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn   = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());

  if( ele->SCluster()->AbsEta() < 1.5 && ( deltaPhiIn > 0.039 || deltaEtaIn > 0.005 ) ) return false;
  else if ( ele->SCluster()->AbsEta() > 1.5 && ( deltaPhiIn > 0.028 || deltaEtaIn > 0.007 ) ) return false;

  return true;
}

Bool_t ElectronTools::PassHggLeptonTagID2012(const Electron* ele) {
  
  if (TMath::Abs(1./ele->E()-1./ele->Pt())>0.05) return false;

  int numInnerHits = ele->Trk()->NExpectedHitsInner();
  if( numInnerHits > 1 ) return false;

  float coviEtaiEta = ele->CoviEtaiEta();
  if( ele->SCluster()->AbsEta() < 1.5 && coviEtaiEta > 0.01 ) return false;
  else if( ele->SCluster()->AbsEta() > 1.5 && coviEtaiEta > 0.03 ) return false; // h

  Double_t deltaPhiIn   = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn   = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());

  if( ele->SCluster()->AbsEta() < 1.5 && ( deltaPhiIn > 0.15 || deltaEtaIn > 0.007 || ele->HadronicOverEm()>0.12) ) return false;   // h
  else if ( ele->SCluster()->AbsEta() > 1.5 && ( deltaPhiIn > 0.10 || deltaEtaIn > 0.009 || ele->HadronicOverEm()>0.10 ) ) return false;   // h

  return true;
}

std::pair<Double_t,Double_t> ElectronTools::ComputeEPCombination( const Electron * ele, const float regression_energy, 
						  const float regression_energy_error) {

  enum Classification { GSF_ELECTRON_UNKNOWN=-1, 
			GSF_ELECTRON_GOLDEN=0, 
			GSF_ELECTRON_BIGBREM=1, 
			GSF_ELECTRON_BADTRACK=2, 
			GSF_ELECTRON_SHOWERING=3, 
			GSF_ELECTRON_GAP=4 } ;

  float newEnergyError_ = regression_energy_error;
  float scEnergy = regression_energy;

  int elClass = ele->Classification();

  float trackMomentum  = ele->PIn() ;
  float errorTrackMomentum_ = 999. ;
  
  // the electron's track momentum error was not available in bambu versions less than 029
  if(ele->TrackMomentumError() > 0)
    errorTrackMomentum_ = ele->TrackMomentumError();
  else if ( ele->GsfTrk()->PtErr() > 0)
    errorTrackMomentum_ = ele->GsfTrk()->PtErr()*cosh(ele->GsfTrk()->Eta());
  else
    assert(0);

  float finalMomentum = ele->E(); // initial
  float finalMomentumError = 999.;
  
  // first check for large errors

  if (errorTrackMomentum_/trackMomentum > 0.5 && regression_energy_error/regression_energy <= 0.5) {
    finalMomentum = regression_energy;    finalMomentumError = regression_energy_error;
  }
  else if (errorTrackMomentum_/trackMomentum <= 0.5 && regression_energy_error/regression_energy > 0.5){  
    finalMomentum = trackMomentum;  finalMomentumError = errorTrackMomentum_;
  }
  else if (errorTrackMomentum_/trackMomentum > 0.5 && regression_energy_error/regression_energy > 0.5){
    if (errorTrackMomentum_/trackMomentum < regression_energy_error/regression_energy) {
      finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;
    }
    else{
      finalMomentum = regression_energy; finalMomentumError = regression_energy_error;
    }
  }
  // then apply the combination algorithm
  else {
     // calculate E/p and corresponding error
    float eOverP = regression_energy / trackMomentum;
    float errorEOverP = sqrt(
			     (regression_energy_error/trackMomentum)*(regression_energy_error/trackMomentum) +
			     (regression_energy*errorTrackMomentum_/trackMomentum/trackMomentum)*
			     (regression_energy*errorTrackMomentum_/trackMomentum/trackMomentum));


    bool eleIsNotInCombination = false ;
    if ( (eOverP  > 1 + 2.5*errorEOverP) || (eOverP  < 1 - 2.5*errorEOverP) || (eOverP < 0.8) || (eOverP > 1.3) )
      { eleIsNotInCombination = true ; }
     if (eleIsNotInCombination)
       {
	 if (eOverP > 1)
	   { finalMomentum = regression_energy ; finalMomentumError = regression_energy_error ; }
	 else
	   {
	     if (elClass == GSF_ELECTRON_GOLDEN)
	       { finalMomentum = regression_energy; finalMomentumError = regression_energy_error; }
	     if (elClass == GSF_ELECTRON_BIGBREM)
	       {
		 if (regression_energy<36)
		   { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		 else
		   { finalMomentum = regression_energy ; finalMomentumError = regression_energy_error ; }
	       }
	     if (elClass == GSF_ELECTRON_BADTRACK)
	       { finalMomentum = regression_energy; finalMomentumError = regression_energy_error ; }
	     if (elClass == GSF_ELECTRON_SHOWERING)
	       {
		 if (regression_energy<30)
		   { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_; }
		 else
		   { finalMomentum = regression_energy; finalMomentumError = regression_energy_error;}
	       }
	     if (elClass == GSF_ELECTRON_GAP)
	       {
		 if (regression_energy<60)
		   { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		 else
		   { finalMomentum = regression_energy; finalMomentumError = regression_energy_error ; }
	       }
	   }
       }     
     else 
       {
	 // combination
	 finalMomentum = (regression_energy/regression_energy_error/regression_energy_error + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
	   (1/regression_energy_error/regression_energy_error + 1/errorTrackMomentum_/errorTrackMomentum_);
	 float finalMomentumVariance = 1 / (1/regression_energy_error/regression_energy_error + 1/errorTrackMomentum_/errorTrackMomentum_);
	 finalMomentumError = sqrt(finalMomentumVariance);
       }
  }  

  std::pair<Double_t,Double_t> value;
  value.first = finalMomentum;
  value.second = finalMomentumError;
  return value;

}
