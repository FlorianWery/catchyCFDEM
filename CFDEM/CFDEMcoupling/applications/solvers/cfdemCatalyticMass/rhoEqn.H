{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(voidfraction,rho)
      + fvc::div(phi)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
    
    Info << "voidfraction max/min: " << max(voidfraction).value() << "/" << min(voidfraction).value() << endl;

}

// ************************************************************************* //
