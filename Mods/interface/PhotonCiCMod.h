//--------------------------------------------------------------------------------------------------
// $Id: PhotonCiCMod.h,v 1.4 2011/07/03 11:04:50 fabstoec Exp $
//
// PhotonCiCMod
//
// This module applies photon identification criteria and exports a pointer to a collection
// of "good photons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONCICMOD_H
#define MITPHYSICS_MODS_PHOTONCICMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

class TNtuple;
class TRandom3;

namespace mithep
{
  class PhotonCiCMod : public BaseMod
  {
    public:
      PhotonCiCMod(const char *name="PhotonCiCMod",
                  const char *title="Photon identification module");

      ~PhotonCiCMod();

      Bool_t              GetApplySpikeRemoval()      const { return fApplySpikeRemoval;   }
      const char         *GetGoodName()               const { return GetGoodPhotonsName(); }
      const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }
      const char         *GetInputName()              const { return fPhotonBranchName;    }
      const char         *GetOutputName()             const { return GetGoodPhotonsName(); }
      Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
      Double_t            GetAbsEtaMax()              const { return fAbsEtaMax;           }
      void                SetApplySpikeRemoval(Bool_t b)    { fApplySpikeRemoval  = b;     }
      void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }
      void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }
      void                SetInputName(const char *n)       { fPhotonBranchName= n;        }
      void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
      void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }
      void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
      void                SetAbsEtaMax(Double_t x)          { fAbsEtaMax       = x;        }
      void                SetPVName(TString s)              { fPVName = s; fPVFromBranch = false; }
      void                SetIsData(Bool_t b)               { fIsData = b; }

      void                AddEnCorrPerRun(UInt_t sRun, UInt_t eRun,
					  Double_t corr_EB_hR9, Double_t corr_EB_lR9,
					  Double_t corr_EE_hR9, Double_t corr_EE_lR9);

      void                SetMCSmearFactors(Double_t _EB_hR9, Double_t _EB_lR9,
					    Double_t _EE_hR9, Double_t _EE_lR9);
      void                FindHiggsPtAndZ(Float_t&, Float_t&);

    protected:
      void                Process();
      void                SlaveBegin();
      unsigned int        FindBestVertex(Photon* ph1, Photon* ph2, const BaseVertex* bsp, bool print=false);
			     
      TString                 fPhotonBranchName;     //name of photon collection (input)
      TString                 fGoodPhotonsName;      //name of exported "good photon" collection
      TString                 fTrackBranchName;      // name of the track collection (only needed for PU corrected isolation)
      TString                 fPileUpDenName;        //name of the PU density collection
      TString                 fElectronName;
      Double_t                fPhotonPtMin;          //min pt cut
      Bool_t                  fApplySpikeRemoval;    //whether apply spike removal
      Double_t                fAbsEtaMax;            //max Abs Eta
      const PhotonCol        *fPhotons;              //!photon branch
      const TrackCol         *fTracks;               //!track branch
      const PileupEnergyDensityCol *fPileUpDen;  //!rho branch
      const ElectronCol      *fElectrons;            //!electron branch

      TString                 fPVName;
      const VertexCol        *fPV;
      Bool_t                  fPVFromBranch;

      const DecayParticleCol *fConversions;
      TString                 fConversionName;

      const BeamSpotCol      *fBeamspot;

      std::vector<Double_t>   fDataEnCorr_EB_hR9;
      std::vector<Double_t>   fDataEnCorr_EB_lR9;
      std::vector<Double_t>   fDataEnCorr_EE_hR9;
      std::vector<Double_t>   fDataEnCorr_EE_lR9;

      std::vector<UInt_t>     fRunStart;
      std::vector<UInt_t>     fRunEnd;

      Double_t                fMCSmear_EB_hR9;
      Double_t                fMCSmear_EB_lR9;
      Double_t                fMCSmear_EE_hR9;
      Double_t                fMCSmear_EE_lR9;

      Bool_t                  fIsData;

      TNtuple*                hCiCTuple;

      TRandom3*               fRnd3;

      TString                 fMCParticleName;
      const MCParticleCol    *fMCParticles;

      TString                 fPileUpName;
      const PileupInfoCol    *fPileUp;

    ClassDef(PhotonCiCMod, 1) // Photon identification module
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::PhotonCiCMod::AddEnCorrPerRun(UInt_t sRun, UInt_t eRun,
						  Double_t corr_EB_hR9, Double_t corr_EB_lR9,
						  Double_t corr_EE_hR9, Double_t corr_EE_lR9)
{
  fDataEnCorr_EB_hR9.push_back(corr_EB_hR9);
  fDataEnCorr_EB_lR9.push_back(corr_EB_lR9);
  fDataEnCorr_EE_hR9.push_back(corr_EE_hR9);
  fDataEnCorr_EE_lR9.push_back(corr_EE_lR9);
  
  fRunStart.push_back(sRun);
  fRunEnd.push_back(eRun);
}

//--------------------------------------------------------------------------------------------------
inline void mithep::PhotonCiCMod::SetMCSmearFactors(Double_t _EB_hR9, Double_t _EB_lR9,
						    Double_t _EE_hR9, Double_t _EE_lR9)
{
  fMCSmear_EB_hR9 = _EB_hR9;
  fMCSmear_EB_lR9 = _EB_lR9;
  fMCSmear_EE_hR9 = _EE_hR9;
  fMCSmear_EE_lR9 = _EE_lR9;
}
#endif
