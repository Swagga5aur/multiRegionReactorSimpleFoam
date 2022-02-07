    fvMesh& mesh = reactorRegions[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields = fieldsReactor[i];
 
    CombustionModel<psiReactionThermo>& reaction = reactionReactor[i];
    psiReactionThermo& thermo = pThermoReactor[i]; 
    PtrList<volScalarField>& Y = reactionReactor[i].thermo().composition().Y(); 
    thermo.validate(args.executable(), "h", "e");
    
    volScalarField& rho = rhoReactor[i];
    volVectorField& U = UReactor[i];
    const dimensionedVector& g = gReactor[i];
    surfaceScalarField& phi = phiReactor[i];

    compressible::turbulenceModel& turb = turbulenceReactor[i];

    volScalarField& K = KReactor[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& Qdot = QdotReactor[i];

    volScalarField& Tc=TcReactor[i];
    volScalarField& zF=zFReactor[i];
    volScalarField& Re=ReReactor[i];
    volScalarField& Nu=NuReactor[i];
    volScalarField& hc=hcReactor[i];
    volScalarField& Pr=PrReactor[i];

    const dimensionedScalar& dp = dpReactor[i]; 
    const dimensionedScalar& eps = epsReactor[i];
    const dimensionedScalar& Kc = KcReactor[i];
    const dimensionedScalar& AoV = AoVReactor[i];

    const dimensionedScalar& a1 = a1Reactor[i];
    const dimensionedScalar& a2 = a2Reactor[i];
    const dimensionedScalar& a3 = a3Reactor[i];

    word cellSetName = "catalyst";

    label zoneID = mesh.cellZones().findZoneID(cellSetName);

    if (zoneID == -1)
    {
	FatalErrorIn("yourFunctionName")
		<< "Cannot find cellZone " << cellSetName << endl
		<< "Valid cellZones are " << mesh.cellZones().names()
		<< exit(FatalError);
    }

    const labelList& cells = mesh.cellZones()[zoneID];

    Info << "Cells in cellzone " << cellSetName << ":" << endl;
    forAll(cells, i)
    {
    const label cell = cells[i];
	zF[cell]=1;
    }

    radiation::radiationModel& rad = radiationReactor[i];

    IOMRFZoneList& MRF = MRFreactor[i];
    //OF4 Update
    fv::options& fvOptions = reactorFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassReactor[i]
    );

    basicMultiComponentMixture& composition = reactionReactor[i].thermo().composition();	

    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << composition.species() << exit(FatalIOError);
    }

    const label inertIndex(composition.species()[inertSpecie]);  
    const bool frozenFlow = frozenFlowReactor[i];
