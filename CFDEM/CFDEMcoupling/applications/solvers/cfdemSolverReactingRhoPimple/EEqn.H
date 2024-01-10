{
    volScalarField& he = thermo.he();
    particleCloud.energyContributions(Qsource);
    particleCloud.energyCoefficients(QCoeff);

    // correct source term (to enable implicit source in EEqn)
    Qsource += QCoeff*(T-he/thermo.Cpv());

    thCond = particleCloud.thermCondM().thermCond();

    fvScalarMatrix EEqn
    (
        fvm::ddt(voidfraction*rho, he) + fvm::div(phi, he)

      + fvc::ddt(voidfraction*rho, K) + fvc::div(phi, K)

      + (
            he.name() == "e"
          ? fvc::div
            (
              fvc::absolute(phi/fvc::interpolate(rho), voidfraction*U),
              p,
              "div(phiv,p)"
            )
          : -dpdt
        )

      - fvm::laplacian
        (
            voidfraction*thCond/thermo.Cpv(),
            he
        )
     ==
        Qdot
      + Qsource + fvm::Sp(QCoeff/thermo.Cpv(), he)
      + fvOptions(voidfraction, rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "Qdot max/min: " << gMax(Qdot) << "/" << gMin(Qdot) << endl;
    Info<< "T max/min: " << gMax(T) << "/" << gMin(T) << endl;

    particleCloud.clockM().start(31,"energySolve");
    particleCloud.solve();
    particleCloud.clockM().stop("energySolve");
}
