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

#ifndef FastJetMassDropTag_h
#define FastJetMassDropTag_h

/** \class FastJetMassDropTag
 *
 *  Finds Boost Jets using mass drop tag algorithm using FastJet library.
 *
 *  \author Yongcheng Wu
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class TIterator;

namespace fastjet
{
class JetDefinition;
class MassDropTagger;

} // namespace fastjet

class FastJetMassDropTag : public DelphesModule
{
public:
  FastJetMassDropTag();
  ~FastJetMassDropTag();

  void Init();
  void Process();
  void Finish();

private:
  fastjet::JetDefinition *fDefinition; //!
  fastjet::MassDropTagger *fMassDropTagger; //!

  Double_t fParameterR;

  Double_t fMassDropTagMu;
  Double_t fMassDropTagYcut;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fConstituentsOutputArray; //!
  TObjArray *fRestConstituentsOutputArray; //!

  ClassDef(FastJetMassDropTag, 1)
};

#endif //FastJetMassDropTag_h
