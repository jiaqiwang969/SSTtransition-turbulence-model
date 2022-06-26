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

Class
    Foam::LESModels::sstTransitionIDDES

Description
    sstTransition IDDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2008).
        A hybrid RANS-LES approach with delayed-DES and wall-modelled LES
        capabilities.
        International Journal of Heat and Fluid Flow, 29(6), 1638-1649.
    \endverbatim

SourceFiles
    sstTransitionIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef sstTransitionIDDES_H
#define sstTransitionIDDES_H

#include "sstTransitionDES.H"
#include "IDDESDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class sstTransitionIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class sstTransitionIDDES
:
    public sstTransitionDES<BasicMomentumTransportModel>
{
    // Private Data

        // Model constants

            dimensionedScalar fwStar_;
            dimensionedScalar cl_;
            dimensionedScalar ct_;

        // Fields

            const IDDESDelta& IDDESDelta_;


    // Private Member Functions

        tmp<volScalarField::Internal> IDDESalpha() const;

        tmp<volScalarField::Internal> ft
        (
            const volScalarField::Internal& magGradU
        ) const;

        tmp<volScalarField::Internal> fl
        (
            const volScalarField::Internal& magGradU
        ) const;

        tmp<volScalarField::Internal> rd
        (
            const volScalarField::Internal& nur,
            const volScalarField::Internal& magGradU
        ) const;

        //- Delay function
        tmp<volScalarField::Internal> fdt
        (
            const volScalarField::Internal& magGradU
        ) const;


protected:

        // Model coefficients

            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;


        using IOobject::modelName;

    
      // Protected Member Functions
        //- Length scale
        virtual tmp<volScalarField::Internal> dTilda
        (
            const volScalarField::Internal& magGradU,
            const volScalarField::Internal& CDES
        ) const;



public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("sstTransitionIDDES");


    // Constructors

        //- Construct from components
        sstTransitionIDDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        sstTransitionIDDES(const sstTransitionIDDES&) = delete;


    //- Destructor
    virtual ~sstTransitionIDDES()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sstTransitionIDDES&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sstTransitionIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
