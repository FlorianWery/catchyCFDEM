if ((!mesh.steady() && !pimple.simpleRho()) || pimple.consistent())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

const volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(voidfraction*rho*rAU));

if (modelType=="A")
{
    rhorAUf *= fvc::interpolate(voidfraction);
}

tmp<volScalarField> rAtU
(
    pimple.consistent()
  ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
  : tmp<volScalarField>(nullptr)
);
tmp<surfaceScalarField> rhorAtUf
(
    pimple.consistent()
  ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(voidfraction*rho*rAtU()))
  : tmp<surfaceScalarField>(nullptr)
);

const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
const surfaceScalarField& rhorAAtUf =
    pimple.consistent() ? rhorAtUf() : rhorAUf;

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiUs("phiUs", fvc::interpolate(voidfraction*rho*rAU*Ksl*Us)& mesh.Sf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(voidfraction*rho)*fvc::flux(HbyA)
);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAAtUf);

if (pimple.transonic())
{
    NotImplemented;
}
else
{
    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(p);
    }
    
    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(voidfraction, rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA)
     ==
        fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + phiUs + pEqn.flux();
        }
    }
}

const bool limitedp = pressureControl.limit(p);

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (limitedp)
{
    rho = thermo.rho();
}

#include "rhoEqn.H"
#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (modelType=="A")
{
    U = HbyA - rAU*(voidfraction*fvc::grad(p)-Ksl*Us);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*Us);
}
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction, p);
}
