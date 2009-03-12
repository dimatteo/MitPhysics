// $Id: JetIDMod.cc,v 1.13 2009/03/03 21:47:34 bendavid Exp $

#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetIDMod)

//--------------------------------------------------------------------------------------------------
JetIDMod::JetIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fJetsName(ModNames::gkPubJetsName),
  fGoodJetsName(ModNames::gkGoodJetsName),  
  fUseJetCorrection(kTRUE),
  fJetPtCut(35.0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void JetIDMod::Process()
{
  // Process entries of the tree. 

  const JetCol  *inJets       = GetObjThisEvt<JetCol>(fJetsName);

  JetOArr *GoodJets = new JetOArr; 
  GoodJets->SetName(fGoodJetsName);

  // loop over jets
  for (UInt_t i=0; i<inJets->GetEntries(); ++i) {
    const Jet *jet = inJets->At(i);

    if (jet->AbsEta() > 5.0) 
      continue;
    
    Double_t jetpt;
    if (fUseJetCorrection)
      jetpt = jet->Pt();
    else
      jetpt = jet->RawMom().Pt();

    if (jetpt < fJetPtCut)
      continue;
    
    // add good jet to collection
    GoodJets->Add(jet);             
  }

  // sort according to pt
  GoodJets->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(GoodJets);  
}

//--------------------------------------------------------------------------------------------------
void JetIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.

}
