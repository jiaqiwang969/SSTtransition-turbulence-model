/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gammaReThetatSSTBase.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
// Tolerance and maximum iteration number for calculation of ReThetat
const scalar gammaReThetatSST::tol_ = 1.0e-4;
const int gammaReThetatSST::maxIter_ = 100;

// Selectable correlation names
const word gammaReThetatSST::CORRN_MENTER2009 = "LangtryMenter2009";
const word gammaReThetatSST::CORRN_SULUKSNA2009 = "SuluksnaEtAl2009";
const word gammaReThetatSST::CORRN_MALAN2009 = "MalanEtAl2009";
const word gammaReThetatSST::CORRN_SORENSEN2009 = "Sorensen2009";
const word gammaReThetatSST::CORRN_TOMAC2013 = "TomacEtAl2013";


// Selectable correlation IDs
const int gammaReThetatSST::CORR_MENTER2009 = 0;
const int gammaReThetatSST::CORR_SULUKSNA2009 = 1;
const int gammaReThetatSST::CORR_MALAN2009 = 2;
const int gammaReThetatSST::CORR_SORENSEN2009 = 3;
const int gammaReThetatSST::CORR_TOMAC2013 = 4;






/*
template<class BasicEddyViscosityModel>
tmp<volScalarField> gammaReThetatSSTBase<BasicEddyViscosityModel>::findGrad
(
    const volVectorField& Vel
) const
{
    volVectorField velgrad(fvc::grad(mag(Vel)));
    
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "dUds",
            Vel.component(0)/(max(mag(Vel),minvel))*velgrad.component(0)
            +Vel.component(1)/max(mag(Vel),minvel)*velgrad.component(1)
        )
    );
}
*/



volScalarField gammaReThetatSST::Flength() const
{
    volScalarField Flength
        (
            IOobject
            (
                "Flength",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            ReThetatTilda_
    );

    switch (corrID_) {
        case CORR_SULUKSNA2009:
            // CORRELATION: SULUKSNA et al. 2009
            Flength = min(
                scalar(0.1)*exp(scalar(-0.022)*ReThetatTilda_+scalar(12))+scalar(0.45),
                scalar(300)
            );
            break;
        case CORR_MALAN2009:
            // CORRELATION: MALAN et al. 2009
            Flength = min(
                exp(scalar(-0.01173)*ReThetatTilda_+scalar(7.168))+scalar(0.5),
                scalar(300)
            );
            break;
        case CORR_SORENSEN2009:
            // CORRELATION: SORENSEN 2009
            Flength = min(
                scalar(150)*exp(scalar(-1)*pow(ReThetatTilda_/scalar(120),1.2))+scalar(0.1),
                scalar(30)
            );
            break;
        case CORR_TOMAC2013:
            // CORRELATION: TOMAC et al. 2013
            Flength = scalar(0.162)+scalar(93.3)*exp(scalar(-1)*sqr(ReThetatTilda_)/scalar(49153))+
                (scalar(50)/(scalar(260)*sqrt(scalar(6.283))))*exp(scalar(-0.5)*
                                                                   sqr((ReThetatTilda_-scalar(520))/scalar(260)));
            break;
        default:
            // CORRELATION: LANGTRY and MENTER 2009
            forAll(Flength, cellI)
            {
                if(ReThetatTilda_[cellI] < scalar(400))
                    Flength[cellI] = scalar(398.189e-1)-scalar(119.270e-4)*ReThetatTilda_[cellI]-scalar(132.567e-6)*sqr(ReThetatTilda_[cellI]);
                else if(ReThetatTilda_[cellI] < scalar(596))
                    Flength[cellI] = scalar(263.404)-scalar(123.939e-2)*ReThetatTilda_[cellI]+scalar(194.548e-5)*sqr(ReThetatTilda_[cellI])-scalar(101.695e-8)*pow3(ReThetatTilda_[cellI]);
                else if(ReThetatTilda_[cellI] < scalar(1200))
                    Flength[cellI] = scalar(0.5)-(ReThetatTilda_[cellI]-scalar(596.0))*scalar(3e-4);
                else
                    Flength[cellI] = scalar(0.3188);
            }

            Flength = (Flength*(scalar(1.0)-exp(-sqr(sqr(y_)*omega_/(0.4*500.0*nu()))))+40.0*exp(-sqr(sqr(y_)*omega_/(0.4*500.0*nu()))));
    }

    return Flength;
}


volScalarField gammaReThetatSST::ReThetac() const
{
    volScalarField ReThetac
        (
            IOobject
            (
                "ReThetac",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            ReThetatTilda_
    );


    switch (corrID_) {
        case CORR_SULUKSNA2009:
            // CORRELATION: SULUKSNA et al. 2009
            ReThetac = min(
                max(
                scalar(1.47)*ReThetatTilda_-sqr(scalar(0.025)*ReThetatTilda_)-scalar(120),
                scalar(125)
                ),
                ReThetatTilda_
            );
            break;
        case CORR_MALAN2009:
            // CORRELATION: MALAN et al. 2009
            ReThetac = min(
                scalar(0.615)*ReThetatTilda_+scalar(61.5),
                ReThetatTilda_
            );
            break;
        case CORR_SORENSEN2009:
            // CORRELATION: SORENSEN 2009
            ReThetac = tanh(pow4((ReThetatTilda_-scalar(100))/scalar(400)))*
                   (ReThetatTilda_+scalar(12000))/scalar(25)+
                   (scalar(1)-tanh(pow4((ReThetatTilda_-scalar(100))/scalar(400))))*
                   (scalar(7)*ReThetatTilda_+scalar(100))/scalar(10);
            break;
        case CORR_TOMAC2013:
            // CORRELATION: TOMAC et al. 2013
            ReThetac = min(
                scalar(0.993)*ReThetatTilda_,
                scalar(0.322)*ReThetatTilda_+(scalar(105900)/(scalar(150)*sqrt(scalar(6.283))))*
                    (exp(scalar(-0.5)*sqr((ReThetatTilda_-scalar(560))/scalar(150)))+
                     exp(scalar(-0.5)*sqr((ReThetatTilda_-scalar(168))/scalar(150))))
            );
            break;
        default:
            // CORRELATION: LANGTRY and MENTER 2009
            forAll(ReThetac, cellI)
            {
                if(ReThetatTilda_[cellI] > scalar(1870))
                    ReThetac[cellI] = ReThetatTilda_[cellI]-(scalar(593.11)+(ReThetatTilda_[cellI]-scalar(1870.0))*scalar(0.482));
                else
                    ReThetac[cellI] = ReThetatTilda_[cellI]-(scalar(396.035e-2)-scalar(120.656e-4)*ReThetatTilda_[cellI]+scalar(868.230e-6)*sqr(ReThetatTilda_[cellI])-scalar(696.506e-9)*pow3(ReThetatTilda_[cellI])+scalar(174.105e-12)*pow4(ReThetatTilda_[cellI]));
            }
    }

    return ReThetac;

}

tmp<volScalarField> gammaReThetatSST::Fonset() const
{
    return tmp<volScalarField>
    (
        new volScalarField
            (
                IOobject
                (
                "Fonset",
                    runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
        max
        (
            min(max(Fonset1(),pow4(Fonset1())),scalar(2))-max(scalar(1)-pow3(Rt()/scalar(2.5)),scalar(0)),
            scalar(0)
        )
        )
    );
}

tmp<volScalarField> gammaReThetatSST::Fturb() const
{
    return exp(-pow4(Rt()/scalar(4)));
}

tmp<volScalarField> gammaReThetatSST::Freattach() const
{
    return exp(-pow4(Rt()/scalar(20)));
}

tmp<volScalarField> gammaReThetatSST::Fwake() const
{
    return exp(-sqr(sqr(y_)*omega_/(scalar(1.0e5)*nu())));
}

tmp<volScalarField> gammaReThetatSST::FThetat() const
{
    volScalarField magVort(sqrt(scalar(2))*mag(skew(fvc::grad(U_))));
    magVort = max(magVort,
                  dimensionedScalar("smallOmega",magVort.dimensions(),SMALL));
    return min
    (
        max
        (
            Fwake()*exp(-pow4(magSqr(U_)
                              /(scalar(375.0)*nu()*magVort*ReThetatTilda_))),
            scalar(1.0)-sqr((ce2_*gamma_-scalar(1.0))/(ce2_-scalar(1.0)))
        ),
        scalar(1.0)
    );
}

// FULL IMPLEMENTATION according to LANGTRY and MENTER 2009
// SULUKSNA et. al. 2009 suggest fixing lambda = 0 and thus omitting
// pressure gradient influence on ReThetat (to be tested!)
void gammaReThetatSST::ReThetat(volScalarField& ReThetatField) const
{
    scalar Tu, lambda, ReThetatOld, ReThetatNew, ReThetatTol, dUds, K;
    volScalarField U2gradU(sqr(U_)&&(fvc::grad(U_)));

    forAll(ReThetatField, cellI)
    {
        int iter = 0;
        Tu = max(
            scalar(100)*sqrt(k_[cellI]/scalar(1.5))/max(mag(U_[cellI]),SMALL),
            scalar(0.027)
        );

        dUds = U2gradU[cellI]/(sqr(max(mag(U_[cellI]),SMALL)));

        // Starting value
        ReThetatNew = max(ReThetatEq(Tu, scalar(0), scalar(0)),scalar(20.0));
        ReThetatTol = ReThetatNew*tol_;

        if (dUds_) {
            do
            {
                ReThetatOld = ReThetatNew;
                lambda = max(
                        min(
                        sqr(ReThetatOld)*nu()()[cellI]*dUds/(sqr(max(mag(U_[cellI]),SMALL))),
                    scalar(0.1)
                    ),
                    scalar(-0.1)
                );
                K = max(
                        min(
                        nu()()[cellI]*dUds/(sqr(max(mag(U_[cellI]),SMALL))),
                    scalar(3e-6)
                    ),
                    scalar(-3e-6)
                );
                ReThetatNew = max(ReThetatEq(Tu, lambda, K),scalar(20.0));

                if (iter++ > maxIter_)
                {
                    FatalErrorIn
                    (
                         "gammaReThetatSST::ReThetat(volScalarField& ReThetatField) const"
                    )   << "Maximum number of iterations exceeded"
                        << abort(FatalError);
                }
            } while(mag(ReThetatNew-ReThetatOld) > ReThetatTol);
        }

        ReThetatField[cellI] = ReThetatNew;
    }

}


scalar gammaReThetatSST::ReThetatEq(scalar Tu, scalar lambda, scalar K) const
{
    scalar FTu;
    scalar FlamK;
    switch (corrID_) {
        case CORR_SULUKSNA2009:
        case CORR_SORENSEN2009:
            // "OLD" CORRELATION FROM MENTER ET AL. (2004)
            FTu = scalar(803.73)*pow((Tu+scalar(0.6067)),scalar(-1.027));
            if(lambda > scalar(0)) {
                scalar FK = scalar(0.0962e6)*K+scalar(0.148e12)*sqr(K)+scalar(0.0141e18)*pow3(K);
                FlamK = scalar(1.0)+FK*(scalar(1.0)-exp(-Tu/scalar(1.5)))+scalar(0.556)*(scalar(1.0)-exp(-scalar(23.9)*lambda))*exp(-Tu/scalar(1.5));
            }
            else {
                scalar Flam = scalar(10.32)*lambda+scalar(89.47)*sqr(lambda)+265.51*pow3(lambda);
                FlamK = scalar(1.0)+Flam*exp(-Tu/scalar(3.0));
            }
            break;
        default:
            // "NEW" CORRELATION FROM LANGTRY/MENTER (2009)
            if(Tu > scalar(1.3)) {
                FTu = scalar(331.5)*pow((Tu-scalar(0.5658)),scalar(-0.671));
            }
            else {
                FTu = scalar(1173.51)-scalar(589.428)*Tu+scalar(0.2196)/sqr(Tu);
            }
            if(lambda > scalar(0)) {
                FlamK = (scalar(1.0)+scalar(0.275)*(scalar(1.0)-exp(scalar(-35.0)*lambda))*exp(scalar(-2.0)*Tu));
            }
            else {
                FlamK = (scalar(1.0)+(scalar(12.986)*lambda+scalar(123.66)*sqr(lambda)+scalar(405.689)*pow3(lambda))*exp(-pow((Tu/scalar(1.5)),scalar(1.5))));
            }
    }
    return FTu*FlamK;
}

tmp<volScalarField> gammaReThetatSST::gammaSep() const
{
    return FThetat()*min
    (
        s1_*Freattach()*max
        (
            sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(scalar(3.235)*nu()*ReThetac())-scalar(1.0),
        scalar(0.0)
        ),
        scalar(2.0)
    );
}

tmp<volScalarField> gammaReThetatSST::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus
    (
        max
        (
            CDkOmega,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        )
    );

    volScalarField arg1
    (
        min
        (
            min
            (
                max
                (
                    (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                    scalar(500)*nu()/(sqr(y_)*omega_)
                ),
                (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
            ),
            scalar(10)
        )
    );

    // Modified blending function!
    return
    max(
        tanh(pow4(arg1)),
    exp(-sqr(pow4(y_*sqrt(k_)/(scalar(120)*nu()))))
    );
}


tmp<volScalarField> gammaReThetatSST::F2() const
{
    volScalarField arg2
    (
        min
        (
            max
            (
                (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            scalar(100)
        )
    );

    return tanh(sqr(arg2));
}

	



/*
template<class BasicEddyViscosityModel>
tmp<volScalarField> gammaReThetatSSTBase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    tmp<volScalarField> R_y = y_ * sqrt(k_) * this->rho_ / this->mu();
    tmp<volScalarField> F3_y = Foam::exp(-Foam::pow(R_y / 120.0, 8.0));
    return max(tanh(pow4(arg1)), F3_y);
    //return tanh(pow4(arg1));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> gammaReThetatSSTBase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> gammaReThetatSSTBase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> gammaReThetatSSTBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void gammaReThetatSSTBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}

template<class BasicEddyViscosityModel> 
void gammaReThetatSSTBase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> gammaReThetatSSTBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
gammaReThetatSSTBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> gammaReThetatSSTBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
	return min
		(
		 GbyNu0,
		 (c1_/a1_)*betaStar_*omega_()
		 *max(a1_*omega_(), b1_*F2*sqrt(S2))
		);
}




template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> gammaReThetatSSTBase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> gammaReThetatSSTBase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> gammaReThetatSSTBase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

*/
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
gammaReThetatSSTBase<BasicEddyViscosityModel>::gammaReThetatSSTBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    ca1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ca1",
            this->coeffDict_,
            2.0
        )
    ),
    ce1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ce1",
            this->coeffDict_,
            1.0
        )
    ),
    ca2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50.0
        )
    ),
    cThetat_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cThetat",
            this->coeffDict_,
            0.03
        )
    ),
    sigmaf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaf",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaThetat_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaThetat",
            this->coeffDict_,
            2.0
        )
    ),
    s1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "s1",
            this->coeffDict_,
            2.0
        )
    ),
    dUds_
    (
        Switch::getOrAddToDict
        (
            "dUds",
            this->coeffDict_,
            false
        )
    ),
    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    kInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kInf",
            coeffDict_,
            0.0,
            sqr(dimLength/dimTime)
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaInf",
            coeffDict_,
            0.0,
            dimless/dimTime
        )
    ),

    y_(wallDist::New(this->mesh_).y()),


    gamma_
    (
        IOobject
        (
            IOobject::groupName("gamma", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    ReThetatTilda_
    (
        IOobject
        (
            IOobject::groupName("ReThetatTilda", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),


    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    nut_
    (
        IOobject
        (
	    IOobject::groupName("nut", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", this->mesh_)
    ),

    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
	(
         "kInf",
         this->coeffDict_,
         k_.dimensions(),
         0
        )
    ),
    omegaInf_
    (
     dimensioned<scalar>::getOrAddToDict
     (
      "omegaInf",
      this->coeffDict_,
      omega_.dimensions(),
      0
     )
    ),
{

    // get correlations name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word corrName
    (
        IOdictionary
        (
            IOobject
            (
                "RASProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookupOrDefault("gammaReThetatSSTCorrelations",CORRN_MENTER2009)
    );


    // set selected correlation ID
    word corrInfo;
    if (corrName == CORRN_SULUKSNA2009) {
        corrID_ = CORR_SULUKSNA2009;
        corrInfo = "Suluksna et al. (2009)";
    }
    else if (corrName == CORRN_MALAN2009) {
        corrID_ = CORR_MALAN2009;
        corrInfo = "Malan et al. (2009)";
    }
    else if (corrName == CORRN_SORENSEN2009) {
        corrID_ = CORR_SORENSEN2009;
        corrInfo = "Sorensen (2009)";
    }
    else if (corrName == CORRN_TOMAC2013) {
        corrID_ = CORR_TOMAC2013;
        corrInfo = "Tomac et al. (2013)";
    }
    else {
        corrID_ = CORR_MENTER2009;
        corrInfo = "Langtry and Menter (2009)";
    }


    Info << "Using gammaReThetat-correlations by " << corrInfo << endl;

    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(scalar(2))*mag(symm(fvc::grad(U_))));
    nut_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class BasicEddyViscosityModel>
void gammaReThetatSSTBase<BasicEddyViscosityModel>::setDecayControl
(
 const dictionary& dict
)
{
	decayControl_.readIfPresent("decayControl", dict);

	if (decayControl_)
	{
		kInf_.read(dict);
		omegaInf_.read(dict);

		Info<< "    Employing decay control with kInf:" << kInf_
			<< " and omegaInf:" << omegaInf_ << endl;
	}
	else
	{
		kInf_.value() = 0;
		omegaInf_.value() = 0;
	}
}
*/

scalar gammaReThetatSST::ReThetatTildaInlet(scalar Tu) const
{
    return ReThetatEq(Tu,scalar(0),scalar(0));
}



tmp<volSymmTensorField> gammaReThetatSST::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> gammaReThetatSST::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}



tmp<fvVectorMatrix> gammaReThetatSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> gammaReThetatSST::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}



template<class BasicEddyViscosityModel>
bool gammaReThetatSSTBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
	ca1_.readIfPresent(coeffDict());
        ce1_.readIfPresent(coeffDict());
        ca2_.readIfPresent(coeffDict());
        ce2_.readIfPresent(coeffDict());
        cThetat_.readIfPresent(coeffDict());
        sigmaf_.readIfPresent(coeffDict());
        sigmaThetat_.readIfPresent(coeffDict());
        s1_.readIfPresent(coeffDict());
        dUds_.readIfPresent("dUds",coeffDict());

        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
	//setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
void gammaReThetatSSTBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField S2(magSqr(symm(fvc::grad(U_))));
    volScalarField G(GName(), nut_*2*S2);
    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();
    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );
    volScalarField F1(this->F1(CDkOmega));


    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)
       *min(2*S2, (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, F2()*sqrt(scalar(2)*S2)))
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
      + beta(F1)*sqr(omegaInf_)
    );

    omegaEqn().relax();







    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    volScalarField nu(this->nu());

    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField str(sqrt(S2));
    volScalarField vort(sqrt(2 * magSqr(skew(tgradU()))));
    volScalarField::Internal GbyNu0
    (
     	this->type() + ":GbyNu",
	(tgradU() && dev(twoSymm(tgradU())))
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu0);
    tgradU.clear();
    


    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

	GbyNu0 = GbyNu(GbyNu0, F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
	  + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
	  + omegaSource()
	  + fvOptions(alpha, rho, omega_)
	);

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)*gamma_
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(min(max(gamma_, 0.1), 1.0) * alpha * rho * betaStar_ * omega_, k_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, k_)
    );
    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);
    correctNut(S2);
//   transition model
    volScalarField Rew(omega_ * y_ * y_ / nu);
    volScalarField Rev(y_ * y_ * str / nu);
    volScalarField Rt(k_ / (omega_ * nu));
    volScalarField F_wake(Foam::exp(-Foam::pow(Rew / 1e5, 2.0)));
    volScalarField delta(50.0 * vort * y_ / max(mag(U), minvel) * 15.0 / 2.0 * nu * ReThetatTilda_ / max(mag(U), minvel));
    volScalarField Rterm1(F_wake * Foam::exp(-Foam::pow(y_ / max(delta, mindis), 4.0)));
    volScalarField Rterm2(1.0 - Foam::pow((gamma_ - 1.0 / 50.0) / (1.0 - 1.0 / 50.0), 2.0));
    volScalarField F_theta(min(max(Rterm1, Rterm2), 1.0));
    volScalarField Ptheta(0.03 * Foam::pow(mag(U), 2.0) / (500.0 * nu) * (1.0 - F_theta, 1.0));
    volScalarField F_turb(Foam::exp(-Foam::pow(Rt / 4., 4.0)));
    volScalarField F_reattach(Foam::exp(-Foam::pow(Rt / 20., 4.0)));
    volScalarField Duds(findGrad(U)); // Velocity gradient along a streamline
    scalar Umag; // absolute value of the velocity
    scalar Tu;   // Turbulent intensity
    rootFunction tf(1, 2, 3, 4);  // creating the function object
    forAll(k_, cellI)
    {
        Umag = max(mag(U[cellI]), 1.e-8); // avoiding division by zero
        Duds[cellI] = max(min(Duds[cellI], 0.5), -0.5); // bounding DUds for robustness
        Tu = 100.0 * sqrt(2.0 / 3.0 * k_[cellI]) / Umag;
        tf.modify(Tu, Duds[cellI], nu[cellI], freeStreamU.value()); // modifying the object
        Reo_[cellI] = NewtonRoot<rootFunction>(tf, 1e-5).root(0.0); // solving the nonlinear equation
        Reo_[cellI] = Reo_[cellI] * freeStreamU.value() / nu[cellI];
    }
    Reo_ = max(Reo_, 20.0); //limiting Reo_ for robustness
 Info << "I am here03"<< endl;
    // Solve the Re_theta equation
    tmp<fvScalarMatrix> ReEqn
    (
        fvm::ddt(alpha, rho, ReThetatTilda_)
      + fvm::div(alphaRhoPhi, ReThetatTilda_)
      - fvm::Sp(fvc::div(alphaRhoPhi), ReThetatTilda_)
      - fvm::laplacian(DRetEff()*alpha*rho, ReThetatTilda_)
     ==
        Ptheta * Reo_*alpha*rho
      - fvm::Sp(Ptheta*alpha*rho, ReThetatTilda_)
    );

    ReEqn.ref().relax();
    solve(ReEqn);

    // Correlation formulas 
    volScalarField Re_crit(ReThetatTilda_);
    volScalarField F_length(ReThetatTilda_);
    forAll(ReThetatTilda_, cellI)
    {
        if (ReThetatTilda_[cellI] <= 1870.0)
        {
            Re_crit[cellI] = ReThetatTilda_[cellI] - 396.035e-02 + 120.656e-04 * ReThetatTilda_[cellI]
                                         - 868.230e-06 * Foam::pow(ReThetatTilda_[cellI], 2.0)
                                         + 696.506e-09 * Foam::pow(ReThetatTilda_[cellI], 3.0)
                                         - 174.105e-12 * Foam::pow(ReThetatTilda_[cellI], 4.0);
        }
        else
        {
            Re_crit[cellI] = ReThetatTilda_[cellI] - 593.11 - 0.482 * (ReThetatTilda_[cellI] - 1870.0);
        }
        if (ReThetatTilda_[cellI] < 400.0)
        {
            F_length[cellI] = 398.189e-01 - 119.270e-04 * ReThetatTilda_[cellI]
                            - 132.567e-06 * Foam::pow(ReThetatTilda_[cellI], 2.0);
        }
        else if (ReThetatTilda_[cellI] >= 400.0 && ReThetatTilda_[cellI] < 596.0)
        {
            F_length[cellI] = 263.404 - 123.939e-02 * ReThetatTilda_[cellI]
                            + 194.548e-05 * Foam::pow(ReThetatTilda_[cellI], 2.0)
                            - 101.695e-08 * Foam::pow(ReThetatTilda_[cellI], 3.0);
        }
         else if (ReThetatTilda_[cellI] >= 596.0 && ReThetatTilda_[cellI] < 1200.0)
        {
            F_length[cellI] = 0.5 - (ReThetatTilda_[cellI] - 596.0) * 3.0e-04;
        }
        else
        {
            F_length[cellI] = 0.3188;
        }
    }
    
    // Solving the intermittency equation
    volScalarField separation_im(min(2.0 * max(0.0, Rev / (3.235 * Re_crit) - 1.0) * F_reattach, 2.0) * F_theta);
    volScalarField F1onset(Rev / (2.193 * max(Re_crit, 1.e-6))); //control the location of the transition onset, change from 2.193 to 3.5 
    volScalarField F2onset(min(max(F1onset, Foam::pow(F1onset, 4.0)), 2.0));
    volScalarField F3onset(max(1.0 - Foam::pow(Rt / 2.5, 3.0), 0.0));
    volScalarField F_onset(max(F2onset - F3onset, 0.0));
    volScalarField F_sub(Foam::exp(-Foam::pow(Rew / (500 * 0.4), 2.0)));
    
    F_length = F_length * (1.0 - F_sub) + 40.0 * F_sub;
    volScalarField Pr1(2.0 * F_length * sqrt(gamma_ * F_onset) * str);
    volScalarField Pr2(0.06 * vort * F_turb * gamma_);
    
    // intermittency equation
    tmp<fvScalarMatrix> imEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::Sp(fvc::div(alphaRhoPhi), gamma_)
      - fvm::laplacian(DimEff()*alpha*rho, gamma_)
    ==
      Pr1*alpha*rho + Pr2*alpha*rho
      - fvm::Sp(Pr1*alpha*rho, gamma_)
      - fvm::Sp(Pr2*alpha*rho * 50.0, gamma_)
    );
    imEqn.ref().relax();
    solve(imEqn);
    gamma_ = max(gamma_, separation_im); //correction for separation induced transition
    gamma_ = min(max(gamma_, 0.00001), 1.0); //bounding intermittency
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
