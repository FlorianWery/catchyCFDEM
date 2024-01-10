// Solve the Momentum equation
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(voidfraction*rho, U) + fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  + fvm::Sp(Ksl, U)
 ==
    voidfraction*rho*g
  + fvOptions(voidfraction, rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn == -fvc::grad(p) + Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
else if (pimple.momentumPredictor())
{
    solve(UEqn == -voidfraction*fvc::grad(p) + Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}

forAll(Us,cellI){
    if(mag(Us[cellI]) == 0)
    {
        PartPresent[cellI] = 0;
    }
    else
    {
        PartPresent[cellI] = 1;
    }
}
