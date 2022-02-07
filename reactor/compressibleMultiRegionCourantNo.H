    scalar CoNum = -GREAT;
    forAll(reactorRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                reactorRegions[regionI],
                runTime,
                rhoReactor[regionI],
                phiReactor[regionI]
            ),
            CoNum
        );
    }
