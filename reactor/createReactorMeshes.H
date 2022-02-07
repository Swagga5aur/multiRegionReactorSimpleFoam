    const wordList reactorNames(rp["reactor"]);

    PtrList<fvMesh> reactorRegions(reactorNames.size());

    forAll(reactorNames, i)
    {
        Info<< "Create reactor mesh for region " << reactorNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        reactorRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    reactorNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
