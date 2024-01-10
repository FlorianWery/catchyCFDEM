Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");
volScalarField& p = thermo.p();

Info<< "Reading field rho" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voidfraction" << endl;
volScalarField voidfraction
(
    IOobject
    (
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field partVolFraction" << endl;
volScalarField partVolFraction
(
    IOobject
    (
        "partVolFraction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0-voidfraction
);

Info<< "Creating fluid-particle heat flux field" << endl;
volScalarField Qsource
(
    IOobject
    (
        "Qsource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
);

Info<< "Creating fluid-particle heat flux coefficient field" << endl;
volScalarField QCoeff
(
    IOobject
    (
        "QCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
);

Info<< "Reading/calculating face flux field phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(rho*U*voidfraction) & mesh.Sf()
);

Info<< "Creating turbulence model" << endl;
autoPtr<fluidThermoPhaseCompressibleMomentumTransportModel> turbulence
(
    fluidThermoPhaseCompressibleMomentumTransportModel::New
    (
        voidfraction,
        rho,
        U,
        phi,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidPhaseThermophysicalTransportModel> thermophysicalTransport
(
    fluidPhaseThermophysicalTransportModel::New(turbulence(), thermo)
);

pressureControl pressureControl
(
    p,
    rho,
    pimple.dict(),
    thermo.incompressible()
);

mesh.setFluxRequired(p.name());

Info<< "Creating field dpdt" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K" << endl;
volScalarField K("K", 0.5*magSqr(U));

Info<< "Reading momentum exchange field Ksl" << endl;
volScalarField Ksl
(
    IOobject
    (
        "Ksl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading particle velocity field Us\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Upartfrac" << endl;
volVectorField Upartfrac
(
    IOobject
    (
        "Upartfrac",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Us*(1.0-voidfraction)
);

Info<< "Reading field PartPresent" << endl;
volScalarField PartPresent
(
    IOobject
    (
        "PartPresent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,0,0,0,0,0,0), 0.0)
);

Info<< "Creating thermal conductivity field thCond" << endl;
volScalarField thCond
(
    IOobject
    (
        "thCond",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1,1,-3,-1,0,0,0), 0.0),
    "zeroGradient"
);
