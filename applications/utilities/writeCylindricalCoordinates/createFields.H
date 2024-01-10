Info << "Reading cell centre coordinates" << endl;
volScalarField ccz
(
	IOobject
	(
	    "ccz",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh.C().component(axDir)
);
volScalarField ccr
(
	IOobject
	(
	    "radial_coordinate",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	sqrt(sqr(mesh.C().component(0))+sqr(mesh.C().component(1)))
);
ccr.write();
volScalarField theta
(
	IOobject
	(
	    "tangential_coordinate",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	atan2(mesh.C().component(1),mesh.C().component(0))
);
theta.write();

if (writeU)
{
    volVectorField U1
    (
        IOobject
        (
            Uname,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField radU
    (
        IOobject
        (
            Uname+".rad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U1.component(1)*mesh.C().component(1)+U1.component(0)*mesh.C().component(0))/ccr
    );
    radU.write();
    volScalarField tanU
    (
        IOobject
        (
            Uname+".tan",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       (U1.component(1)*mesh.C().component(0)-U1.component(0)*mesh.C().component(1))/ccr
    );
    tanU.write();
}

if (writeUPrime)
{
    volSymmTensorField U1
    (
        IOobject
        (
            UPrimename,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField radUVar
    (
        IOobject
        (
            UPrimename+".rad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U1.component(0)*sqr(mesh.C().component(0))+U1.component(3)*sqr(mesh.C().component(1))
            +2*U1.component(1)*mesh.C().component(0)*mesh.C().component(1))/sqr(ccr)
    );
    radUVar.write();
    volScalarField tanUVar
    (
        IOobject
        (
            UPrimename+".tan",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U1.component(0)*sqr(mesh.C().component(1))+U1.component(3)*sqr(mesh.C().component(0))
            -2*U1.component(1)*mesh.C().component(0)*mesh.C().component(1))/sqr(ccr)
    );
    tanUVar.write();
}