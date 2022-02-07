// Initialise reactor field pointer lists
PtrList<psiReactionThermo> pThermoReactor(reactorRegions.size());
PtrList<CombustionModel<psiReactionThermo>> reactionReactor(reactorRegions.size());
PtrList<volScalarField> rhoReactor(reactorRegions.size());
PtrList<volScalarField> TcReactor(reactorRegions.size());
PtrList<volScalarField> zFReactor(reactorRegions.size());
PtrList<volVectorField> UReactor(reactorRegions.size());
PtrList<surfaceScalarField> phiReactor(reactorRegions.size());
PtrList<uniformDimensionedVectorField> gReactor(reactorRegions.size());
PtrList<volScalarField> CpReactor(reactorRegions.size());
PtrList<compressible::turbulenceModel> turbulenceReactor(reactorRegions.size());
PtrList<radiation::radiationModel> radiationReactor(reactorRegions.size());
PtrList<volScalarField> KReactor(reactorRegions.size());
PtrList<volScalarField> QdotReactor(reactorRegions.size());

List<scalar> initialMassReactor(reactorRegions.size());
List<bool> frozenFlowReactor(reactorRegions.size(), false);
PtrList<dimensionedScalar> rhoMaxReactor(reactorRegions.size());
PtrList<dimensionedScalar> rhoMinReactor(reactorRegions.size());

PtrList<dimensionedScalar> epsReactor(reactorRegions.size());
PtrList<dimensionedScalar> KcReactor(reactorRegions.size());
PtrList<dimensionedScalar> dpReactor(reactorRegions.size());
PtrList<dimensionedScalar> AoVReactor(reactorRegions.size());

PtrList<dimensionedScalar> a1Reactor(reactorRegions.size());
PtrList<dimensionedScalar> a2Reactor(reactorRegions.size());
PtrList<dimensionedScalar> a3Reactor(reactorRegions.size());

PtrList<volScalarField> ReReactor(reactorRegions.size());
PtrList<volScalarField> PrReactor(reactorRegions.size());
PtrList<volScalarField> NuReactor(reactorRegions.size());
PtrList<volScalarField> hcReactor(reactorRegions.size());

//OF4 Update
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsReactor(reactorRegions.size());

PtrList<IOMRFZoneList> MRFreactor(reactorRegions.size());
PtrList<fv::options> reactorFvOptions(reactorRegions.size());

// Populate reactor field pointer lists
forAll(reactorRegions, i)
{
    Info<< "*** Reading reactor mesh thermophysical properties for region "
        << reactorRegions[i].name() << nl << endl;

    Info<< "Creating reaction model\n" << endl;

    pThermoReactor.set
    (
        i,
          autoPtr<psiReactionThermo>(psiReactionThermo::New(reactorRegions[i])).ptr()
    );

    Info<< "    Adding to rhoReactor\n" << endl;
    rhoReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pThermoReactor[i].rho()
        )
    );

    Info<< "    Adding to UReactor\n" << endl;
    UReactor.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i]
        )
    );

    Info<< "    Adding to phiReactor\n" << endl;
    phiReactor.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoReactor[i]*UReactor[i])
                & reactorRegions[i].Sf()
        )
    );

    Info<< "    Adding to gReactor\n" << endl;
    gReactor.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                reactorRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    CpReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Cp",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pThermoReactor[i].Cp()

        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulenceReactor.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoReactor[i],
            UReactor[i],
            phiReactor[i],
            pThermoReactor[i]
        ).ptr()
    );

    Info<< "Creating reaction model\n" << endl;
 
    reactionReactor.set
    (
        i,
          autoPtr<CombustionModel<psiReactionThermo>>
          (
          CombustionModel<psiReactionThermo>::New
           (
            pThermoReactor[i],
            turbulenceReactor[i]          
           )
          )
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            reactorRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading porousity eps\n" << endl;
    epsReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("eps")
	)
    );

    Info<< "Reading catalyst conductivity\n" << endl;
    KcReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("Kc")
	)
    );

    Info<< "Reading pellet diameter\n" << endl;
    dpReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("dp")
	)
    );

    Info<< "Reading pellet diameter\n" << endl;
    a1Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a1")
	)
    );

    Info<< "Reading pellet diameter\n" << endl;
    a2Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a2")
	)
    );

    Info<< "Reading pellet diameter\n" << endl;
    a3Reactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("a3")
	)
    );

    Info<< "Reading pellet diameter\n" << endl;
    AoVReactor.set
    (
	i,
        new dimensionedScalar
        (
       		transportProperties.lookup("AoV")
	)
    );


    Info<< "    Adding to ReReactor\n" << endl;
    ReReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Re",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i],
	    dimensionedScalar("Re", dimless, scalar(0.0))
//	    rhoReactor[i]*mag(UReactor[i])*dpReactor[i]/turbulenceReactor[i].muEff()
        )
    );

    Info<< "    Adding to PrReactor\n" << endl;
    PrReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Pr",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i],
	    dimensionedScalar("Pr", dimless, scalar(0.7))
//	    rhoReactor[i]*mag(UReactor[i])*dpReactor[i]/turbulenceReactor[i].muEff()
        )
    );

    Info<< "    Adding to NuReactor\n" << endl;
    NuReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Nu",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i],
	    dimensionedScalar("Nu", dimless, scalar(0.0))
        )
    );

    Info<< "    Adding to ReReactor\n" << endl;
    hcReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "hc",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i],
	    dimensionedScalar("hc", dimensionSet(1, 0, -3, -1, 0, 0, 0), scalar(0.0))
        )
    );

    Info<< "    Adding to TcReactor\n" << endl;
    TcReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Tc",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i]
        )
    );

    Info<< "    Adding to zFReactor\n" << endl;
    zFReactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "zF",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
	    reactorRegions[i],
	    dimensionedScalar("zF", dimless, scalar(0.0))
        )
    );

    radiationReactor.set
    (
        i,
        radiation::radiationModel::New(reactionReactor[i].thermo().T())
    );

    initialMassReactor[i] = fvc::domainIntegrate(rhoReactor[i]).value();

    Info<< "    Adding to KReactor\n" << endl;
    KReactor.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UReactor[i])
        )
    );

   Info<< "    Adding to QdotReactor\n" << endl;
   QdotReactor.set
   (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                reactorRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
                reactorRegions[i],
                dimensionedScalar
            (
                "Qdot",
                dimEnergy/dimVolume/dimTime,
                0
            )
        )
    );

    const dictionary& simpleDict =
        reactorRegions[i].solutionDict().subDict("SIMPLE");

    simpleDict.readIfPresent("frozenFlow", frozenFlowReactor[i]);

    rhoMaxReactor.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMaxReactor",
                simpleDict,
                dimDensity,
                GREAT
            )
        )
    );

    rhoMinReactor.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMinReactor",
                simpleDict,
                dimDensity,
                0
            )
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFreactor.set
    (
        i,
        new IOMRFZoneList(reactorRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
//OF4 Update
    reactorFvOptions.set
     (
        i,
        new fv::options(reactorRegions[i])
     );
 
    fieldsReactor.set
    	(
    	    i,
			new multivariateSurfaceInterpolationScheme<scalar>::fieldTable 
		);

    basicSpecieMixture& composition = reactionReactor[i].thermo().composition();		
	
	forAll(composition.Y(), j)
	{
		    fieldsReactor[i].add(composition.Y()[i]);
	}


		
    fieldsReactor[i].add(reactionReactor[i].thermo().he());
    turbulenceReactor[i].validate();    
}
