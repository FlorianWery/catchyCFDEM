/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "WenYuDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WenYuDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    WenYuDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
WenYuDrag::WenYuDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    minVoidfraction_(propsDict_.lookupOrDefault<scalar>("minVoidfraction",0.1)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce");    // first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");         // other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          // other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction"); // other are debug
    particleCloud_.probeM().writeHeader();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(0,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(0,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(0,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(0,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WenYuDrag::~WenYuDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WenYuDrag::setForce() const
{
    if (scaleDia_ > 1)
    {
        Info << "WenYuDrag using scale = " << scaleDia_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "WenYuDrag using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
    label cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            dragExplicit = vector::zero;
            voidfraction=0;
            dragCoefficient=0;
            Ufluid = vector::zero;

            if (cellI > -1) // particle Found
            {
                if(forceSubM(0).interpolation())
                {
                    position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;

                    if (voidfraction > 1.00) voidfraction = 1.0;
                    if (voidfraction < minVoidfraction_) voidfraction = minVoidfraction_;
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                ds = 2*particleCloud_.radius(index);
                nuf = nufField[cellI];
                rho = rhoField[cellI];
                magUr = mag(Ur);
                Rep = 0;
              

                if (magUr > 0)
                {

                    // calc particle Re Nr and 
                    Rep = ds/scaleDia_*voidfraction*magUr/(nuf+SMALL);

                    // calc particle's drag coefficient (i.e., Force per unit slip velocity and Stokes drag)
                    dragCoefficient = F(voidfraction, Rep)
                                      *3.0*M_PI*nuf*rho*voidfraction*ds
                                      *pow(voidfraction, -3.65)
                                      *scaleDia_
                                      *scaleDrag_;
                    if (modelType_=="B")
                        dragCoefficient /= voidfraction;

                    drag = dragCoefficient*Ur; //total drag force!

                    forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],forceSubM(0).verbose(),index);
                }

                if(forceSubM(0).verbose() && index >-1 && index <102)
                {
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "drag (total) = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
        }
}

double WenYuDrag::F(double voidfraction, double Rep) const
{

    double F0 = 0.0;

    if (Rep < 1000)  
    {
      F0 = (1.0 + 0.15*Foam::pow(Rep, 0.687));
    }
    else
    {
      F0 = (0.44/24.0)*Rep;
    }
    return   F0;

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
