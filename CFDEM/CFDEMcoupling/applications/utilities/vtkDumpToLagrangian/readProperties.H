// Create properties file object
IOdictionary propsDict
(
    IOobject
    (
        "vtkDumpProperties",
        runTime.system(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

word dump_filename(fileName(propsDict.lookup("demDumpFiles")).expand());
word selector(fileName(propsDict.lookup("selector")));
scalar demTimeStep(readScalar(propsDict.lookup("demTimeStep")));

scalar demTimeStepNumber(readScalar(propsDict.lookup("demTimeStepNumber")));

word cloudName_(propsDict.lookupOrDefault<word>("cloudName","particleCloud"));

List<word> scalarFlds;
if (propsDict.found("scalarFields"))
{
    propsDict.lookup("scalarFields")() >> scalarFlds;
}

List<word> vectorFlds;
if (propsDict.found("vectorFields"))
{
    propsDict.lookup("vectorFields")() >> vectorFlds;
}
