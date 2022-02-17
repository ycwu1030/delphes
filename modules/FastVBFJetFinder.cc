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

/** \class FastVBFJetFinder
 *
 *  Finds VBF jets using FastJet library.
 *
 *  \author Yongcheng Wu
 *
 */

#include "modules/FastVBFJetFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//------------------------------------------------------------------------------

FastVBFJetFinder::FastVBFJetFinder() :
  fPlugin(0), fRecomb(0), fAxesDef(0), fMeasureDef(0), fNjettinessPlugin(0), fValenciaPlugin(0),
  fDefinition(0), fAreaDefinition(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

FastVBFJetFinder::~FastVBFJetFinder()
{
}

//------------------------------------------------------------------------------

void FastVBFJetFinder::Init()
{
  JetDefinition::Plugin *plugin = 0;
  JetDefinition::Recombiner *recomb = 0;

  // define algorithm

  fJetAlgorithm = GetInt("JetAlgorithm", 6);
  fParameterR = GetDouble("ParameterR", 0.5);

  fConeRadius = GetDouble("ConeRadius", 0.5);
  fSeedThreshold = GetDouble("SeedThreshold", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fMaxIterations = GetInt("MaxIterations", 100);
  fMaxPairSize = GetInt("MaxPairSize", 2);
  fIratch = GetInt("Iratch", 1);
  fAdjacencyCut = GetInt("AdjacencyCut", 2);
  fOverlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetPTMin = GetDouble("JetPTMin", 30.0);
  fJetEtaMin = GetDouble("JetEtaMin", 2.0);
  // fJetPairDeltaEtaMin = GetDouble("JetPairDeltaEtaMin", 5.0);
  // fJetPairInvariantMassMin = GetDouble("JetPairInvariantMassMin", 500)

  switch(fJetAlgorithm)
  {
  case 1:
    plugin = new CDFJetCluPlugin(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
    fDefinition = new JetDefinition(plugin);
    break;
  case 2:
    plugin = new CDFMidPointPlugin(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
    fDefinition = new JetDefinition(plugin);
    break;
  case 3:
    plugin = new SISConePlugin(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
    fDefinition = new JetDefinition(plugin);
    break;
  case 4:
    fDefinition = new JetDefinition(kt_algorithm, fParameterR);
    break;
  case 5:
    fDefinition = new JetDefinition(cambridge_algorithm, fParameterR);
    break;
  default:
  case 6:
    fDefinition = new JetDefinition(antikt_algorithm, fParameterR);
    break;
  case 7:
    recomb = new WinnerTakeAllRecombiner();
    fDefinition = new JetDefinition(antikt_algorithm, fParameterR, recomb, Best);
    break;
  }

  fPlugin = plugin;
  fRecomb = recomb;

  ClusterSequence::print_banner();

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
  fConstituentsOutputArray = ExportArray(GetString("ConstituentsOutputArray", "constituents"));
  fRestConstituentsOutputArray = ExportArray(GetString("RestConstituentsOutputArray", "rests"));
}

//------------------------------------------------------------------------------

void FastVBFJetFinder::Finish()
{
  vector<TEstimatorStruct>::iterator itEstimators;

  for(itEstimators = fEstimators.begin(); itEstimators != fEstimators.end(); ++itEstimators)
  {
    if(itEstimators->estimator) delete itEstimators->estimator;
  }

  if(fItInputArray) delete fItInputArray;
  if(fDefinition) delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fPlugin) delete static_cast<JetDefinition::Plugin *>(fPlugin);
  if(fRecomb) delete static_cast<JetDefinition::Recombiner *>(fRecomb);
  if(fNjettinessPlugin) delete static_cast<JetDefinition::Plugin *>(fNjettinessPlugin);
  if(fAxesDef) delete fAxesDef;
  if(fMeasureDef) delete fMeasureDef;
  if(fValenciaPlugin) delete static_cast<JetDefinition::Plugin *>(fValenciaPlugin);
}

//------------------------------------------------------------------------------

void FastVBFJetFinder::Process()
{
  Candidate *candidate, *constituent;
  TLorentzVector momentum;

  Double_t deta, dphi, detaMax, dphiMax;
  Double_t time, timeWeight;
  Double_t neutralEnergyFraction, chargedEnergyFraction;

  Int_t number, ncharged, nneutrals;
  Int_t charge;
  Double_t rho = 0.0;
  PseudoJet jet, area;
  ClusterSequence *sequence;
  vector<PseudoJet> inputList, outputList, subjets, tmpList;
  vector<PseudoJet>::iterator itInputList, itOutputList, itTmpList1, itTmpList2;
  unordered_set<Int_t> VBFConstituentID;
  vector<TEstimatorStruct>::iterator itEstimators;
  Double_t excl_ymerge23 = 0.0;
  Double_t excl_ymerge34 = 0.0;
  Double_t excl_ymerge45 = 0.0;
  Double_t excl_ymerge56 = 0.0;

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

  sequence = new ClusterSequence(inputList, *fDefinition);

  outputList.clear();

  outputList = sorted_by_pt(sequence->inclusive_jets(fJetPTMin));

  // loop over all jets and see if any of them have eta larger than requirement
  // which will be stored into tmpList and called as forward jet
  tmpList.clear();
  for(itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList)
  {
    jet = *itOutputList;
    if(fJetAlgorithm == 7) jet = join(jet.constituents());

    if(abs(jet.eta()) > fJetEtaMin)
    {
      tmpList.push_back(jet);
    }
  }
  // Check which two of the possible forward jet will be paired to be VBF jet
  int forward_jet_size = tmpList.size();
  int forward_jet_id1 = -1;
  int forward_jet_id2 = -1;
  detaMax = 0.0;
  for(int id1 = 0; id1 < forward_jet_size; id1++)
  {
    for(int id2 = id1 + 1; id2 < forward_jet_size; id2++)
    {
      double eta1 = tmpList[id1].eta();
      double eta2 = tmpList[id2].eta();
      if(eta1 * eta2 > 0) continue;
      deta = abs(eta1 - eta2);
      if(deta > detaMax)
      {
        detaMax = deta;
        forward_jet_id1 = id1;
        forward_jet_id2 = id2;
      }
    }
  }

  // loop over VBF jets and export them
  detaMax = 0.0;
  dphiMax = 0.0;

  VBFConstituentID.clear();
  for(int fid = 0; fid < forward_jet_size; ++fid)
  {
    if(fid != forward_jet_id1 && fid != forward_jet_id2) continue;
    jet = tmpList[fid];
    if(fJetAlgorithm == 7) jet = join(jet.constituents());

    momentum.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());

    area.reset(0.0, 0.0, 0.0, 0.0);
    if(fAreaDefinition) area = itOutputList->area_4vector();

    candidate = factory->NewCandidate();

    time = 0.0;
    timeWeight = 0.0;

    charge = 0;

    ncharged = 0;
    nneutrals = 0;

    neutralEnergyFraction = 0.;
    chargedEnergyFraction = 0.;

    inputList.clear();
    inputList = sequence->constituents(jet);

    for(itInputList = inputList.begin(); itInputList != inputList.end(); ++itInputList)
    {
      if(itInputList->user_index() < 0) continue;
      constituent = static_cast<Candidate *>(fInputArray->At(itInputList->user_index()));
      VBFConstituentID.insert(itInputList->user_index());

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
    if(VBFConstituentID.count(cid) != 0) continue;
    constituent = static_cast<Candidate *>(fInputArray->At(cid));
    fRestConstituentsOutputArray->Add(constituent);
  }

  delete sequence;
}
