# SSTtransition-turbulence-model
- Aim: Implantation of SSTtransition turbulence model 

- Ref: [SSTtransition turbulence model](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/#YEAR_2020) with [newversion](https://www.cfd-online.com/Forums/openfoam-solving/180356-sst-transition.html)



## [Validate in v2012](https://github.com/jiaqiwang969/Axis-2Dbump):
Improve: Original code is [sstTransition](https://gitlab.com/tilasoldo/openfoam_share/-/blob/master/src/TurbulenceModels/turbulenceModels/Base/sstTransition/sstTransitionBase.C), update to ***compressible*** and combine it to IDDES. Tip: multiply "\*alpha\*rho" in solving equation


### How to use?

Add turbulenceProperties
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


## V9 version
Validate？

<img src="https://cdn.mathpix.com/snip/images/Zg9KfN9V065uPpVOvtGGXDVDeCuHH91U53jO1BYY8mo.original.fullsize.png" width="340px">

