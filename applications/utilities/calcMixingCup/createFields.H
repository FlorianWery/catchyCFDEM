Info << "Reading face centre coordinates" << endl;
surfaceScalarField ccz
(
	IOobject
	(
	    "ccz",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	mesh.Cf().component(axDir)
);

labelList order;
scalarList cZ(ccz.internalField());
scalarList biZ(ccz.boundaryField()[inPatch]);
scalarList boZ(ccz.boundaryField()[outPatch]);
sortedOrder(cZ, order);
scalar outZ=average(boZ);
scalar inZ=average(biZ);
scalar deltaZ =(outZ-inZ)/nPoints;
Info << "Min axial position = " << inZ << ", max axial position = " << outZ << ", deltaZ for averaging = " << deltaZ<< endl;
List<scalarList> iFields(fieldIDs.size());
List<scalarList> biFields(fieldIDs.size());
List<scalarList> boFields(fieldIDs.size());

Info << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
		"U",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    ),
    mesh
);

Info << "Reading field phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
		"phi",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    ),
    mesh
);

forAll(fieldIDs, fi)
{
    Info << "Reading field " << fieldIDs[fi] << endl;
    volScalarField volFieldI
    (
        IOobject
        (
            word(fieldIDs[fi]),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info << "Converting field " << fieldIDs[fi] << " to surfaceScalarField" << endl;
    surfaceScalarField fieldI
    (
        IOobject
        (
            word(fieldIDs[fi]),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(volFieldI)
    );

    Info << "Getting internal field of " << fieldIDs[fi] << endl;
    iFields[fi]=fieldI.internalField();

    Info << "Getting boundary fields of " << fieldIDs[fi] << endl;
    biFields[fi]=fieldI.boundaryField()[inPatch];
    boFields[fi]=fieldI.boundaryField()[outPatch];
}
