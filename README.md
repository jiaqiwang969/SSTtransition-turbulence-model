# SSTtransition-turbulence-model
- Aim: Parametric analysis for  SSTtransition turbulence model
<<<<<<< HEAD
- Ref: [SSTtransition turbulence model](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/#YEAR_2020) with [newversion](https://www.cfd-online.com/Forums/openfoam-solving/180356-sst-transition.html)



## Modified for v2012:
Original code is [sstTransition](https://gitlab.com/tilasoldo/openfoam_share/-/blob/master/src/TurbulenceModels/turbulenceModels/Base/sstTransition/sstTransitionBase.C), update to ***compressible*** and combine it to IDDES. 

Tip: multiply "*alpha*rho" in solving equation


```
simulationType  LES;

LES
{
    LESModel        sstTransitionIDDES;
    printCoeffs     no;
    turbulence      yes;
    delta           IDDESDelta;
    IDDESDeltaCoeffs
    {
        hmax           maxDeltaxyzCubeRoot;
        maxDeltaxyzCubeRootCoeffs
        {
        }
    }
}

```


Todo: for v9
=======
- Ref: [SSTtransition turbulence model](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/#YEAR_2020) with [newversion](https://www.cfd-online.com/Forums/openfoam-solving/180356-sst-transition.html) [gammaReTheta](https://gitlab.com/tilasoldo/openfoam_share/-/tree/master/src)
>>>>>>> 2867806ac48d28e6904eadbc07295aa0a5d62288
