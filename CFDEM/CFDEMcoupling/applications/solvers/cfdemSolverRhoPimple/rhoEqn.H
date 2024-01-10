{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(voidfraction,rho)
      + fvc::div(phi)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
}

// ************************************************************************* //
