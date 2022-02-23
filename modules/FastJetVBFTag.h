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

#ifndef FastJetVBFTag_h
#define FastJetVBFTag_h

/** \class FastJetVBFTag
 *
 *  Finds VBF jets using FastJet library.
 *
 *  \author Yongcheng Wu
 *
 */

#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;
class TIterator;

namespace fastjet
{
class JetDefinition;
class AreaDefinition;
class JetMedianBackgroundEstimator;
class PseudoJet;
namespace contrib
{
class NjettinessPlugin;
class ValenciaPlugin;
class AxesDefinition;
class MeasureDefinition;
} // namespace contrib
} // namespace fastjet

typedef Double_t (*TVBFMethod_ptr)(const fastjet::PseudoJet &, const fastjet::PseudoJet &);

class FastJetVBFTag : public DelphesModule
{
public:
  FastJetVBFTag();
  ~FastJetVBFTag();

  void Init();
  void Process();
  void Finish();

private:
  void *fPlugin; //!
  void *fRecomb; //!

  fastjet::JetDefinition *fDefinition; //!

  Int_t fJetAlgorithm;
  Double_t fParameterR;

  Double_t fJetPTMin;
  Double_t fJetEtaMin;
  Int_t fVBFMethod;

  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;
  Int_t fAdjacencyCut;
  Double_t fOverlapThreshold;

  TVBFMethod_ptr fMethod;
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fConstituentsOutputArray; //!
  TObjArray *fRestConstituentsOutputArray; //!

  ClassDef(FastJetVBFTag, 1)
};

#endif
