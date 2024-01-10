/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PFA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class GThermoType, class SThermoType>
Foam::chemistryGSReductionMethods::PFA<CompType, GThermoType, SThermoType>::PFA
(
    const IOdictionary& dict,
    TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
)
:
    chemistryGSReductionMethod<CompType, GThermoType, SThermoType>(dict, chemistry),
    searchInitGasSet_(this->coeffsDict_.subDict("initialGasSet").size()),
    searchInitSolidSet_(this->coeffsDict_.subDict("initialSolidSet").size()),
    searchInitSet_(this->coeffsDict_.subDict("initialGasSet").size() + this->coeffsDict_.subDict("initialSolidSet").size())
{
  label j=0;
  dictionary initGasSet = this->coeffsDict_.subDict("initialGasSet");
  for (label i=0; i < this->nGasSpecie_; i++)
  {
      if (initGasSet.found(chemistry.Y()[i].member()))
      {
          searchInitGasSet_[j]=i;
          searchInitSet_[j++]=i;
      }
  }
  if (j<searchInitGasSet_.size())
  {
      FatalErrorInFunction
          << searchInitGasSet_.size()-j
          << " species in the initial set is not in the mechanism "
          << initGasSet
          << exit(FatalError);
  }

  j=0;
  dictionary initSolidSet = this->coeffsDict_.subDict("initialSolidSet");
  for (label i=0; i < this->nSolidSpecie_; i++)
  {
      if (initSolidSet.found(chemistry.Ys()[i].member()))
      {
          searchInitSolidSet_[j]=i;
          searchInitSet_[j + searchInitGasSet_.size()]=i + this->nGasSpecie_;
          j++;
      }
  }
  if (j<searchInitSolidSet_.size())
  {
      FatalErrorInFunction
          << searchInitSolidSet_.size()-j
          << " species in the initial set is not in the mechanism "
          << initSolidSet
          << exit(FatalError);
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class GThermoType, class SThermoType>
Foam::chemistryGSReductionMethods::PFA<CompType, GThermoType, SThermoType>::~PFA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class GThermoType, class SThermoType>
void Foam::chemistryGSReductionMethods::PFA<CompType, GThermoType, SThermoType>::reduceMechanism
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
)
{
    scalarField& completeC(this->chemistry_.completeC());
    scalarField c1(this->chemistry_.nEqns(), 0.0);

    for (label i=0; i<this->nSpecie_; i++)
    {
        c1[i] = c[i];
        completeC[i] = c[i];
    }

    c1[this->nSpecie_] = T;
    c1[this->nSpecie_+1] = p;

    // Compute the rAB matrix
    RectangularMatrix<scalar> PAB(this->nSpecie_,this->nSpecie_,0.0);
    RectangularMatrix<scalar> CAB(this->nSpecie_,this->nSpecie_,0.0);
    scalarField PA(this->nSpecie_,0.0);
    scalarField CA(this->nSpecie_,0.0);

    // Number of initialized rAB for each lines
    Field<label> NbrABInit(this->nSpecie_,0);
    // Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos(this->nSpecie_, this->nSpecie_, -1);
    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie_, this->nSpecie_, -1);

    scalar pf, cf, pr, cr;
    label lRef, rRef;
    forAll(this->chemistry_.gasReactions(), i)
    {
        const Reaction<GThermoType>& R = this->chemistry_.gasReactions()[i];
        // for each reaction compute omegai
        scalar omegai = R.omega
        (
            p, T, c1, li, pf, cf, lRef, pr, cr, rRef
        );

        // then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)

        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());

        forAll(R.lhs(), s)// compute rAB for all species in the left hand side
        {
            label ss = R.lhs()[s].index;
            scalar sl = -R.lhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        forAll(R.rhs(), s) // compute rAB for all species in the right hand side
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff; // vAi = v''-v' => here v''
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        wAID.shrink();
        forAll(wAID, id)
        {
            label curID = wAID[id];
            scalar curwA = wA[id];
            List<bool> deltaBi(this->nSpecie_, false);
            FIFOStack<label> usedIndex;
            forAll(R.lhs(),j)
            {
                label sj = R.lhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.rhs(),j)
            {
                label sj = R.rhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }

            deltaBi[curID] = false;

            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();

                if (deltaBi[curIndex])
                {
                    deltaBi[curIndex] = false;
                    if (rABPos(curID, curIndex)==-1)
                    {
                        rABPos(curID, curIndex) = NbrABInit[curID];
                        rABOtherSpec(curID, NbrABInit[curID]) = curIndex;
                        if (curwA > 0.0)
                        {
                            PAB(curID, NbrABInit[curID]) = curwA;
                        }
                        else
                        {
                            CAB(curID, NbrABInit[curID]) = -curwA;
                        }
                        NbrABInit[curID]++;
                    }
                    else
                    {
                        if (curwA > 0.0)
                        {
                            PAB(curID, rABPos(curID, curIndex)) += curwA;
                        }
                        else
                        {
                            CAB(curID, rABPos(curID, curIndex)) += -curwA;
                        }
                    }
                }
            }
            // Now that every species of the reactions has been visited, we can
            // compute the production and consumption rate. It avoids getting
            // wrong results when species are present in both lhs and rhs
            if (curwA > 0.0)
            {
                if (PA[curID] == 0.0)
                {
                    PA[curID] = curwA;
                }
                else
                {
                    PA[curID] += curwA;
                }
            }
            else
            {
                if (CA[curID] == 0.0)
                {
                    CA[curID] = -curwA;
                }
                else
                {
                    CA[curID] += -curwA;
                }
            }
        }
    }

    forAll(this->chemistry_.solidReactions(), i)
    {
        const SurfaceReaction<SThermoType>& R = this->chemistry_.solidReactions()[i];
        // for each reaction compute omegai
        scalar omegai = R.omega
        (
            p, T, c1, li, pf, cf, lRef, pr, cr, rRef
        );

        // then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)

        DynamicList<scalar> wA(R.lhs().size()+R.glhs().size()+R.rhs().size()+R.grhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.glhs().size()+R.rhs().size()+R.grhs().size());

        forAll(R.lhs(), s)// compute rAB for all species in the left hand side
        {
            label ss = R.lhs()[s].index + this->nGasSpecie_;
            scalar sl = -R.lhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        forAll(R.glhs(), s)// compute rAB for all species in the left hand side
        {
            label ss = R.glhs()[s].index;
            scalar sl = -R.glhs()[s].stoichCoeff; // vAi = v''-v' => here -v'
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        forAll(R.rhs(), s) // compute rAB for all species in the right hand side
        {
            label ss = R.rhs()[s].index + this->nGasSpecie_;
            scalar sl = R.rhs()[s].stoichCoeff; // vAi = v''-v' => here v''
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        forAll(R.grhs(), s) // compute rAB for all species in the right hand side
        {
            label ss = R.grhs()[s].index;
            scalar sl = R.grhs()[s].stoichCoeff; // vAi = v''-v' => here v''
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        wAID.shrink();
        forAll(wAID, id)
        {
            label curID = wAID[id];
            scalar curwA = wA[id];
            List<bool> deltaBi(this->nSpecie_, false);
            FIFOStack<label> usedIndex;
            forAll(R.lhs(),j)
            {
                label sj = R.lhs()[j].index + this->nGasSpecie_;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.glhs(),j)
            {
                label sj = R.glhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.rhs(),j)
            {
                label sj = R.rhs()[j].index + this->nGasSpecie_;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }
            forAll(R.grhs(),j)
            {
                label sj = R.grhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }

            deltaBi[curID] = false;

            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();

                if (deltaBi[curIndex])
                {
                    deltaBi[curIndex] = false;
                    if (rABPos(curID, curIndex)==-1)
                    {
                        rABPos(curID, curIndex) = NbrABInit[curID];
                        rABOtherSpec(curID, NbrABInit[curID]) = curIndex;
                        if (curwA > 0.0)
                        {
                            PAB(curID, NbrABInit[curID]) = curwA;
                        }
                        else
                        {
                            CAB(curID, NbrABInit[curID]) = -curwA;
                        }
                        NbrABInit[curID]++;
                    }
                    else
                    {
                        if (curwA > 0.0)
                        {
                            PAB(curID, rABPos(curID, curIndex)) += curwA;
                        }
                        else
                        {
                            CAB(curID, rABPos(curID, curIndex)) += -curwA;
                        }
                    }
                }
            }
            // Now that every species of the reactions has been visited, we can
            // compute the production and consumption rate. It avoids getting
            // wrong results when species are present in both lhs and rhs
            if (curwA > 0.0)
            {
                if (PA[curID] == 0.0)
                {
                    PA[curID] = curwA;
                }
                else
                {
                    PA[curID] += curwA;
                }
            }
            else
            {
                if (CA[curID] == 0.0)
                {
                    CA[curID] = -curwA;
                }
                else
                {
                    CA[curID] += -curwA;
                }
            }
        }
    }
    // rii = 0.0 by definition

    // Compute second generation link strength
    // For all species A look at all rAri of the connected species ri and
    // compute rriB with all the connected species of ri, B different of A.  If
    // a new species is connected, add it to the list of connected species.  It
    // is a connection of second generation and it will be aggregated in the
    // final step to evaluate the total connection strength (or path flux).
    // Compute rsecond=rAri*rriB with A!=ri!=B
    RectangularMatrix<scalar> PAB2nd(this->nSpecie_,this->nSpecie_,0.0);
    RectangularMatrix<scalar> CAB2nd(this->nSpecie_,this->nSpecie_,0.0);

    // Number of initialized rAB for each lines
    Field<label> NbrABInit2nd(this->nSpecie_, 0);

    // Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos2nd(this->nSpecie_, this->nSpecie_, -1);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec2nd
    (
        this->nSpecie_, this->nSpecie_, -1
    );

    forAll(NbrABInit, A)
    {
        for (int i=0; i<NbrABInit[A]; i++)
        {
            label ri = rABOtherSpec(A, i);
            scalar maxPACA = max(PA[ri],CA[ri]);
            if (maxPACA > vSmall)
            {
                for (int j=0; j<NbrABInit[ri]; j++)
                {
                    label B = rABOtherSpec(ri, j);
                    if (B != A) // if B!=A and also !=ri by definition
                    {
                        if (rABPos2nd(A, B)==-1)
                        {
                            rABPos2nd(A, B) = NbrABInit2nd[A]++;
                            rABOtherSpec2nd(A, rABPos2nd(A, B)) = B;
                            PAB2nd(A, rABPos2nd(A, B)) =
                                PAB(A, i)*PAB(ri, j)/maxPACA;
                            CAB2nd(A, rABPos2nd(A, B)) =
                                CAB(A, i)*CAB(ri, j)/maxPACA;
                        }
                        else
                        {
                            PAB2nd(A, rABPos2nd(A, B)) +=
                                PAB(A, i)*PAB(ri, j)/maxPACA;
                            CAB2nd(A, rABPos2nd(A, B)) +=
                                CAB(A, i)*CAB(ri, j)/maxPACA;
                        }
                    }
                }
            }
        }
    }

    // Using the rAB matrix (numerator and denominator separated)
    label speciesNumber = 0;
    label speciesGasNumber = 0;
    label speciesSolidNumber = 0;

    // set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }

    // Initialize the FIFOStack for search set
    const labelList& SIS(this->searchInitSet_);
    FIFOStack<label> Q;

    for (label i=0; i<SIS.size(); i++)
    {
        label q = SIS[i];
        this->activeSpecies_[q] = true;
        speciesNumber++;
        if(q < this->nGasSpecie_)
        {
            speciesGasNumber++;
        }
        else
        {
            speciesSolidNumber++;
        }
        Q.push(q);
    }

    // Execute the main loop for R-value
    while (!Q.empty())
    {
        label u = Q.pop();
        scalar Den = max(PA[u],CA[u]);

        if (Den!=0.0)
        {
            // first generation
            for (label v=0; v<NbrABInit[u]; v++)
            {
                label otherSpec = rABOtherSpec(u, v);
                scalar rAB = (PAB(u, v)+CAB(u, v))/Den; // first generation
                label id2nd = rABPos2nd(u, otherSpec);
                if (id2nd !=-1)// if there is a second generation link
                {
                    rAB += (PAB2nd(u, id2nd)+CAB2nd(u, id2nd))/Den;
                }
                // the link is stronger than the user-defined tolerance
                if
                (
                    rAB >= this->tolerance()
                 && !this->activeSpecies_[otherSpec]
                )
                {
                    Q.push(otherSpec);
                    this->activeSpecies_[otherSpec] = true;
                    speciesNumber++;
                    if(otherSpec < this->nGasSpecie_)
                    {
                        speciesGasNumber++;
                    }
                    else
                    {
                        speciesSolidNumber++;
                    }
                }

            }
            // second generation link only (for those without first link)
            for (label v=0; v<NbrABInit2nd[u]; v++)
            {
                label otherSpec = rABOtherSpec2nd(u, v);
                scalar rAB = (PAB2nd(u, v)+CAB2nd(u, v))/Den;
                // the link is stronger than the user-defined tolerance
                if
                (
                    rAB >= this->tolerance()
                 && !this->activeSpecies_[otherSpec]
                )
                {
                    Q.push(otherSpec);
                    this->activeSpecies_[otherSpec] = true;
                    speciesNumber++;
                    if(otherSpec < this->nGasSpecie_)
                    {
                        speciesGasNumber++;
                    }
                    else
                    {
                        speciesSolidNumber++;
                    }
                }
            }
        }
    }

    // Put a flag on the reactions containing at least one removed species

    forAll(this->chemistry_.gasReactions(), i)
    {
        const Reaction<GThermoType>& R = this->chemistry_.gasReactions()[i];
        this->chemistry_.gasReactionsDisabled()[i] = false;

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;

            // The species is inactive then the reaction is removed
            if (!this->activeSpecies_[ss])
            {
                // Flag the reaction to disable it
                this->chemistry_.gasReactionsDisabled()[i] = true;
                break;
            }
        }

        // If the reaction has not been disabled yet
        if (!this->chemistry_.gasReactionsDisabled()[i])
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!this->activeSpecies_[ss])
                {
                    this->chemistry_.gasReactionsDisabled()[i] = true;
                    break;
                }
            }
        }
    }

    forAll(this->chemistry_.solidReactions(), i)
    {
        const SurfaceReaction<SThermoType>& R = this->chemistry_.solidReactions()[i];
        this->chemistry_.solidReactionsDisabled()[i] = false;

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index + this->nGasSpecie_;

            // The species is inactive then the reaction is removed
            if (!this->activeSpecies_[ss])
            {
                // Flag the reaction to disable it
                this->chemistry_.solidReactionsDisabled()[i] = true;
                break;
            }
        }

        if (!this->chemistry_.solidReactionsDisabled()[i])
        {
            forAll(R.glhs(), s)
            {
                label ss = R.glhs()[s].index;
                if (!this->activeSpecies_[ss])
                {
                    this->chemistry_.solidReactionsDisabled()[i] = true;
                    break;
                }
            }
            if (!this->chemistry_.solidReactionsDisabled()[i])
            {
                forAll(R.rhs(), s)
                {
                    label ss = R.rhs()[s].index + this->nGasSpecie_;
                    if (!this->activeSpecies_[ss])
                    {
                        this->chemistry_.solidReactionsDisabled()[i] = true;
                        break;
                    }
                }

                if (!this->chemistry_.solidReactionsDisabled()[i])
                {
                    forAll(R.grhs(), s)
                    {
                        label ss = R.grhs()[s].index;
                        if (!this->activeSpecies_[ss])
                        {
                            this->chemistry_.solidReactionsDisabled()[i] = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    this->NsSimp_ = speciesNumber;
    this->NsGasSimp_ = speciesGasNumber;
    this->NsSolidSimp_ = speciesSolidNumber;
    scalarField& simplifiedC(this->chemistry_.simplifiedC());
    simplifiedC.setSize(this->NsSimp_+2);
    DynamicList<label>& s2c(this->chemistry_.simplifiedToCompleteIndex());
    s2c.setSize(this->NsSimp_);
    Field<label>& c2s(this->chemistry_.completeToSimplifiedIndex());

    label j = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        if (this->activeSpecies_[i])
        {
            s2c[j] = i;
            simplifiedC[j] = c[i];
            c2s[i] = j++;
            if (!this->chemistry_.active(i))
            {
                this->chemistry_.setActive(i);
            }
        }
        else
        {
            c2s[i] = -1;
        }
    }
    simplifiedC[this->NsSimp_] = T;
    simplifiedC[this->NsSimp_+1] = p;
    this->chemistry_.setNsDAC(this->NsSimp_);
    // change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->NsSimp_);
    this->chemistry_.setNGasSpecie(this->NsGasSimp_);
    this->chemistry_.setNSolidSpecie(this->NsSolidSimp_);
}


// ************************************************************************* //
