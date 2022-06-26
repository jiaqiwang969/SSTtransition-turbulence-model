/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "sstTransitionIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::IDDESalpha() const
{
    return volScalarField::Internal::New
    (
        modelName("alpha"),
        max(0.25 - this->y_()/IDDESDelta_.hmax(), scalar(-5))
    );
}






template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::ft
(
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        modelName("ft"),
        tanh(pow3(sqr(Ct_)*rd(this->nut_, magGradU)))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::fl
(
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        modelName("fl"),
        tanh(pow(sqr(Cl_)*rd(this->nu(), magGradU), 10))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::rd
(
    const volScalarField::Internal& nur,
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        modelName("rd"),
        min
        (
            nur
           /(
               max
               (
                   magGradU,
                   dimensionedScalar(magGradU.dimensions(), small)
               )*sqr(this->kappa_*this->y_())
            ),
            scalar(10)
        )
    );
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::fdt
(
     const volScalarField::Internal& magGradU
) const
{
    return 1 - tanh(pow(Cdt1_*rd(this->nut_, magGradU), Cdt2_));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
sstTransitionIDDES<BasicMomentumTransportModel>::dTilda
(
    const volScalarField::Internal& magGradU,
    const volScalarField::Internal& CDES
) const
{
    const volScalarField::Internal& k = this->k_;
    const volScalarField::Internal& omega = this->omega_;

    const volScalarField::Internal lRAS(sqrt(k)/(this->betaStar_*omega));
    const volScalarField::Internal lLES(CDES*this->delta());

    const volScalarField::Internal alpha(IDDESalpha());
    //const volScalarField alpha(this->alpha());
    //const volScalarField expTerm(exp(sqr(alpha)));
    const volScalarField::Internal expTerm
    (
        modelName("expTerm"),
        exp(sqr(alpha))
    );


    tmp<volScalarField::Internal> fB = min(2*pow(expTerm, -9.0), scalar(1));
    tmp<volScalarField::Internal> fe1 =
        2*(pos0(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));
    tmp<volScalarField::Internal> fe2 = 1 - max(ft(magGradU), fl(magGradU));
    tmp<volScalarField::Internal> fe = max(fe1 - 1, scalar(0))*fe2;

    const volScalarField::Internal fdTilda(max(1 - fdt(magGradU), fB));

/*    return max
    (
        fdTilda*(1 + fe)*lRAS + (1 - fdTilda)*lLES,
        dimensionedScalar(dimLength, small)
    );
*/
    return volScalarField::Internal::New
    (
        modelName("dTilda"),
        max
        (
            dimensionedScalar(dimLength, small),
	    fdTilda*(1 + fe)*lRAS + (1 - fdTilda)*lLES
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
sstTransitionIDDES<BasicMomentumTransportModel>::sstTransitionIDDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    sstTransitionDES<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    Cdt1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdt1",
            this->coeffDict_,
            20
        )
    ),
    Cdt2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdt2",
            this->coeffDict_,
            3
        )
    ),
    Cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            this->coeffDict_,
            5
        )
    ),
    Ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            this->coeffDict_,
            1.87
        )
    ),
    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool sstTransitionIDDES<BasicMomentumTransportModel>::read()
{
    if (sstTransitionDES<BasicMomentumTransportModel>::read())
    {
	Cdt1_.readIfPresent(this->coeffDict());
        Cdt2_.readIfPresent(this->coeffDict());
        Cl_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
