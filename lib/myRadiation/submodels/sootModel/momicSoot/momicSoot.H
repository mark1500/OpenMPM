/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::radiation::momicSoot

Description
    This soot model is purely an state model. The ammount of soot produced is
    determined by a single step chemistry as :

        nuf Fuel + nuOx Ox = nuP P + nuSoot soot

    nuSoot is prescribed by the user.

    The single step chemistry used is read from the combustion.
    The soot is not considered into the thermodynamics of the system and it
    is not considered as an extra specie in the solver.

    The spacial distribution is given by the normalization of the first product
    on the rhs of the reaction by default or it can be added as input.

    The input dictionary reads like in the radiationProperties dictionary:

    sootModel momicSoot<gasHThermoPhysics>;

    momicSootCoeffs
    {
        nuSoot              0.015;
        Wsoot               12;
        mappingField        P;
    }

SourceFiles
    momicSoot.C

\*---------------------------------------------------------------------------*/

#ifndef momicSoot_H
#define momicSoot_H

#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"

#include "fluidThermo.H"
#include "reactingMixture.H"
#include "turbulentFluidThermoModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class momicSoot Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class momicSoot
:
    public sootModel
{

    // Static functions

        //- Check mixture in thermo
        static const reactingMixture<ThermoType>& checkThermo
        (
            const fluidThermo&
        );


    // Private data

      //- Soot M0 (M0_ = real M0/rho_)
      volScalarField M0_;

      //- Soot M1 (M1_ = real M1/rho_)
      volScalarField M1_;

      //-Soot M2 (M2_ = real M2/rho_)
      volScalarField M2_;

      //-Soot M3 (M3_ = real M3/rho_)
      volScalarField M3_;

      //-Soot M4 (M4_ = real M4/rho_)
      volScalarField M4_;

      //- Soot model dictionary
      dictionary coeffsDict_;

      //- Inception rate for real M0
      volScalarField incepRateM0_;

      //- Inception rate for real M1
      volScalarField incepRateM1_;

      //- Inception rate for real M2
      volScalarField incepRateM2_;

      //- Inception rate for real M3
      volScalarField incepRateM3_;

      //- Inception rate for real M4
      volScalarField incepRateM4_;

     	 // constat parameters
      scalar Sct_;
     	scalar rhoS_;
     	scalar NA_;
      scalar cperPAH_;
      scalar cperDimer_;
      scalar diamPAH_;

        //- Thermo package
        const fluidThermo& thermo_;

public:

    //- Runtime type information
    TypeName("momicSoot");


    // Constructors

        //- Construct from components
        momicSoot
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~momicSoot();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();


        // Access

            //- Return Ysoot //is this one necessary?
            const volScalarField& soot() const
            {
                return M0_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "momicSoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
