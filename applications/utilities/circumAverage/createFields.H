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
	    "ccr",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	sqrt(sqr(mesh.C().component((axDir+2)%3))+sqr(mesh.C().component((axDir+1)%3)))
);

ccr.write();
DynamicList<label> labelConstantZ;
DynamicList<scalar> radPosZ;
forAll(ccz, i)
{
    if(ccz[i]>planeZ-deltaZ && ccz[i]<planeZ+deltaZ)
    {
        labelConstantZ.append(i);
        radPosZ.append(ccr[i]);
    }
}

labelList orderR;
sortedOrder(radPosZ, orderR);
scalar maxR = radPosZ[orderR[orderR.size()-1]];
scalar minR = radPosZ[orderR[0]];
scalar deltaR =(maxR-minR)/nPoints;
Info << "Min radial position = " << minR << ", max radial position = " << maxR << ", delta for averaging = " << deltaR<< endl;

scalarList cv(labelConstantZ.size());
forAll(labelConstantZ,ri)
{
    cv[ri]=mesh.V()[labelConstantZ[ri]];
}
    
List<scalarList> iFields(fieldIDs.size());

Info << "Reading field U1" << endl;
volVectorField U1
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("U",phaseName1)+"Mean" :
            IOobject::groupName("U",phaseName1)
        ),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField radU1_
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("radU",phaseName1)+"Mean" :
            IOobject::groupName("radU",phaseName1)
        ),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (U1.component((axDir+2)%3)*mesh.C().component((axDir+2)%3)+U1.component((axDir+1)%3)*mesh.C().component((axDir+1)%3))/ccr
);
radU1_.write();
volScalarField tanU1_
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("tanU",phaseName1)+"Mean" :
            IOobject::groupName("tanU",phaseName1)
        ),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::sqrt(Foam::sqr(mag(U1))-Foam::sqr(radU1_))
);
tanU1_.write();
scalarList magU1(labelConstantZ.size());
scalarList radU1(labelConstantZ.size());
scalarList tanU1(labelConstantZ.size());
forAll(labelConstantZ,ri)
{
    magU1[ri]=mag(U1[labelConstantZ[ri]]);
    radU1[ri]=radU1_[labelConstantZ[ri]];
    tanU1[ri]=tanU1_[labelConstantZ[ri]];
}

/*Info << "Reading field U2" << endl;
volVectorField U2
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("U",phaseName2)+"Mean" :
            IOobject::groupName("U",phaseName2)
        ),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField radU2_
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("radU",phaseName2)+"Mean" :
            IOobject::groupName("radU",phaseName2)
        ),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (U2.component((axDir+2)%3)*mesh.C().component((axDir+2)%3)+U2.component((axDir+1)%3)*mesh.C().component((axDir+1)%3))/ccr
);
radU2_.write();
volScalarField tanU2_
(
    IOobject
    (
        (
            args.optionFound("meanValue")
            ?
            IOobject::groupName("tanU",phaseName2)+"Mean" :
            IOobject::groupName("tanU",phaseName2)
        ),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::sqrt(Foam::sqr(mag(U2))-Foam::sqr(radU2_))
);
tanU2_.write();
scalarList magU2(labelConstantZ.size());
scalarList radU2(labelConstantZ.size());
scalarList tanU2(labelConstantZ.size());
forAll(labelConstantZ,ri)
{
    magU2[ri]=mag(U2[labelConstantZ[ri]]);
    radU2[ri]=radU2_[labelConstantZ[ri]];
    tanU2[ri]=tanU2_[labelConstantZ[ri]];
}*/

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
    
    scalarList fieldI(labelConstantZ.size());
    forAll(labelConstantZ,ri)
    {
        fieldI[ri]=volFieldI[labelConstantZ[ri]];
    }
    iFields[fi]=fieldI;
}
