particleCloud.clockM().start(29,"Y");

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

if (solveSpecies)
{
    reaction->correct();
    Qdot = voidfraction*reaction->Qdot();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(voidfraction, rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)

              - fvm::laplacian
                (
                    voidfraction*turbulence->muEff(),
                    Yi
                )
              ==
                voidfraction*reaction->R(Yi)
              + fvOptions(voidfraction, rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    if (inertIndex != -1)
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
        /*forAll(Y, i)
        {
            if (i != inertIndex)
            {
                Y[i] *= (scalar(1)-Y[inertIndex])/Yt;
            }
        }*/
    }
    else
    {
        forAll(Y, i)
        {
            Y[i] /= Yt;
        }
    }
}
particleCloud.clockM().stop("Y");
