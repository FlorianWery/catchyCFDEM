if (modelType=="A")
{
    volScalarField rUsed = rUA*voidfraction;
    constrainPressure(p, U, phi, rUsed);
}
else constrainPressure(p, U, phi, rUA);
