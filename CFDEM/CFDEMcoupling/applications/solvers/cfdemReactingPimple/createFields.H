Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

label inertIndex = -1;
if (thermo.found("inertSpecie"))
{
    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << composition.species() << exit(FatalIOError);
    }
    inertIndex = composition.species()[inertSpecie];
}

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

autoPtr<rhoReactionPhaseThermophysicalTransportModel>
thermophysicalTransport
(
    rhoReactionPhaseThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

pressureControl pressureControl
(
    p,
    rho,
    pimple.dict(),
    thermo.incompressible()
);

mesh.setFluxRequired(p.name());

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> reaction
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);

Info<< "Creating field dpdt\n" << endl;
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

Info<< "Creating field kinetic energy K\n" << endl;
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

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
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
