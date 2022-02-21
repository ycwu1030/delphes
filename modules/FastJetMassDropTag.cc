/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class FastJetMassDropTag
 *
 *  Finds Boost Jets using mass drop tag algorithm using FastJet library.
 *
 *  \author Yongcheng Wu
 *
 */

#include "modules/FastJetMassDropTag.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TString.h"

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"

using namespace std;
using namespace fastjet;

//------------------------------------------------------------------------------

FastJetMassDropTag::FastJetMassDropTag() :
  fDefinition(0), fItInputArray(0), fMassDropTagger(0) {}

//------------------------------------------------------------------------------

FastJetMassDropTag::~FastJetMassDropTag() {}

//------------------------------------------------------------------------------

void FastJetMassDropTag::Init()
{
  fParameterR = GetDouble("ParameterR", 1.2);

  fMassDropTagMu = GetDouble("MassDropTagMu", 0.667);
  fMassDropTagYcut = GetDouble("MassDropTagYcut", 0.09);

  fDefinition = new JetDefinition(cambridge_algorithm, fParameterR);
  fMassDropTagger = new MassDropTagger(fMassDropTagMu, fMassDropTagYcut);

  ClusterSequence::print_banner();

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  fConstituentsOutputArray = ExportArray(GetString("ConstituentsOutputArray", "constituents"));
  fRestConstituentsOutputArray = ExportArray(GetString("RestConstituentsOutputArray", "rests"));
}

//------------------------------------------------------------------------------

void FastJetMassDropTag::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fDefinition) delete fDefinition;
  if(fMassDropTagger) delete fMassDropTagger;
}

//------------------------------------------------------------------------------

void FastJetMassDropTag::Process()
{

  Candidate *candidate, *constituent;
  TLorentzVector momentum;

  Double_t deta, dphi, detaMax, dphiMax;
  Double_t time, timeWeight;
  Double_t neutralEnergyFraction, chargedEnergyFraction;
  Int_t ncharged, nneutrals;
  Int_t charge;
  Double_t rho = 0.0;
  Double_t excl_ymerge23 = 0.0;
  Double_t excl_ymerge34 = 0.0;
  Double_t excl_ymerge45 = 0.0;
  Double_t excl_ymerge56 = 0.0;

  Int_t number;

  PseudoJet jet, area;
  vector<PseudoJet> inputList;
  unordered_set<Int_t> ConstituentID;

  DelphesFactory *factory = GetFactory();

  inputList.clear();

  // loop over input objects
  fItInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;
  }

  // * fDefinition will always be cambridge_algorithm but with user provided dR
  ClusterSequence *sequence = new ClusterSequence(inputList, *fDefinition);

  vector<PseudoJet> jets = sorted_by_pt(sequence->inclusive_jets());

  // * Now perform the jet tagging using mass drop tagger to the hardest one
  PseudoJet tagged_jet = (*fMassDropTagger)(jets[0]);
  if(tagged_jet == 0)
  {
    // * No substructure found, just return;
    delete sequence;
    return;
  }

  // * When substructure is found, then "filter" it
  vector<PseudoJet> sub_jets = tagged_jet.pieces();
  PseudoJet sub_jet_1 = sub_jets[0];
  PseudoJet sub_jet_2 = sub_jets[1];
  double Rsubs = sub_jet_1.delta_R(sub_jet_2);
  double Rfilt = min(Rsubs / 2.0, 0.3);
  unsigned nfilt = 3;
  Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(nfilt));
  PseudoJet filtered_jet = filter(tagged_jet);

  vector<PseudoJet> outputList = filtered_jet.pieces();
  for(auto itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList)
  {
    jet = *itOutputList;

    momentum.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());
    area.reset(0.0, 0.0, 0.0, 0.0);
    candidate = factory->NewCandidate();

    time = 0.0;
    timeWeight = 0.0;

    charge = 0;

    ncharged = 0;
    nneutrals = 0;

    neutralEnergyFraction = 0.0;
    chargedEnergyFraction = 0.0;

    inputList.clear();
    inputList = sequence->constituents(jet);
    for(auto itInputList = inputList.begin(); itInputList != inputList.end(); ++itInputList)
    {
      if(itInputList->user_index() < 0) continue;
      constituent = static_cast<Candidate *>(fInputArray->At(itInputList->user_index()));
      ConstituentID.insert(itInputList->user_index());

      deta = TMath::Abs(momentum.Eta() - constituent->Momentum.Eta());
      dphi = TMath::Abs(momentum.DeltaPhi(constituent->Momentum));
      if(deta > detaMax) detaMax = deta;
      if(dphi > dphiMax) dphiMax = dphi;

      if(constituent->Charge == 0)
      {
        nneutrals++;
        neutralEnergyFraction += constituent->Momentum.E();
      }
      else
      {
        ncharged++;
        chargedEnergyFraction += constituent->Momentum.E();
      }

      time += TMath::Sqrt(constituent->Momentum.E()) * (constituent->Position.T());
      timeWeight += TMath::Sqrt(constituent->Momentum.E());

      charge += constituent->Charge;

      fConstituentsOutputArray->Add(constituent);
      candidate->AddCandidate(constituent);
    }

    candidate->Momentum = momentum;
    candidate->Position.SetT(time / timeWeight);
    candidate->Area.SetPxPyPzE(area.px(), area.py(), area.pz(), area.E());

    candidate->DeltaEta = detaMax;
    candidate->DeltaPhi = dphiMax;
    candidate->Charge = charge;
    candidate->NNeutrals = nneutrals;
    candidate->NCharged = ncharged;

    candidate->NeutralEnergyFraction = (momentum.E() > 0) ? neutralEnergyFraction / momentum.E() : 0.0;
    candidate->ChargedEnergyFraction = (momentum.E() > 0) ? chargedEnergyFraction / momentum.E() : 0.0;

    //for exclusive clustering, access y_n,n+1 as exclusive_ymerge (fNJets);
    candidate->ExclYmerge23 = excl_ymerge23;
    candidate->ExclYmerge34 = excl_ymerge34;
    candidate->ExclYmerge45 = excl_ymerge45;
    candidate->ExclYmerge56 = excl_ymerge56;

    fOutputArray->Add(candidate);
  }
  for(int cid = 0; cid < fInputArray->GetEntries(); ++cid)
  {
    if(ConstituentID.count(cid) != 0) continue;
    constituent = static_cast<Candidate *>(fInputArray->At(cid));
    fRestConstituentsOutputArray->Add(constituent);
  }

  delete sequence;
}
