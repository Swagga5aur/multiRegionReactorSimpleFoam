{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(eps*turb.alphaEff(), he)
//      - fvm::laplacian(eps*turb.muEff(), he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
//      + reaction.Sh()
      + fvOptions(rho, he)
      - fvc::Sp(hc*AoV,thermo.T())*zF
      + hc*AoV*Tc*zF
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    rad.correct();

//    thermo.T().min(dimensionedScalar("thermo.T()",thermo.T().dimensions(),1000));
//    thermo.T().max(dimensionedScalar("thermo.T()",thermo.T().dimensions(),373));

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
