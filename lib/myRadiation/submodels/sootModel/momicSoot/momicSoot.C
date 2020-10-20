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

\*---------------------------------------------------------------------------*/

#include "momicSoot.H"
#include "reactingMixture.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


template<class ThermoType>
const Foam::reactingMixture<ThermoType>&
Foam::radiation::momicSoot<ThermoType>::checkThermo
(
    const fluidThermo& thermo
)
{
    if (isA<reactingMixture<ThermoType>>(thermo))
    {
        return dynamic_cast<const reactingMixture<ThermoType>& >
        (
            thermo
        );
    }
    else
    {
        FatalErrorInFunction
            << "Inconsistent thermo package for " << thermo.type()
            << "Please select a thermo package based on "
            << "reactingMixture" << exit(FatalError);

        return dynamic_cast<const reactingMixture<ThermoType>& >
        (
            thermo
        );
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::momicSoot<ThermoType>::momicSoot
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),

    M0_
    (
        IOobject
        (
            "M0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    M1_
    (
        IOobject
        (
            "M1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    M2_
    (
        IOobject
        (
            "M2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    M3_
    (
        IOobject
        (
            "M3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    M4_
    (
        IOobject
        (
            "M4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),

    incepRateM0_
    (
      IOobject
      (
        "incepRateM0",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("incepRateM0", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    incepRateM1_
    (
      IOobject
      (
        "incepRateM1",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("incepRateM1", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    incepRateM2_
    (
      IOobject
      (
        "incepRateM2",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("incepRateM2", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    incepRateM3_
    (
      IOobject
      (
        "incepRateM3",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("incepRateM3", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    incepRateM4_
    (
      IOobject
      (
        "incepRateM4",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("incepRateM4", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    // initialize the parameters used in the momic soot model
    Sct_(readScalar(coeffsDict_.lookup("Sct"))),
    rhoS_(readScalar(coeffsDict_.lookup("rhoS"))), //unit: kg/m3
    NA_(readScalar(coeffsDict_.lookup("NA"))),
    cperPAH_(readScalar(coeffsDict_.lookup("cperPAH"))), //unit:#
    cperDimer_(readScalar(coeffsDict_.lookup("cperDimer"))), //unit:#
    diamPAH_(readScalar(coeffsDict_.lookup("diamPAH"))), //unit:m
    //
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::momicSoot<ThermoType>::~momicSoot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiation::momicSoot<ThermoType>::correct()
{
    // define the momic soot model: transportation equations
    // read the flux, temperature and pressure used in the soot model
    const surfaceScalarField& phi_(mesh_.lookupObject<surfaceScalarField>("phi"));
    const volScalarField& rho_(mesh_.lookupObject<volScalarField>("rho")); //unit:kg/m3
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    // define some chemistry constants
    scalar kb = 1.3806488e-23; // BOLTZMANN constant, unit: m^2 kg s^(-2) K^(-1)
    scalar cMass = 2.004e-26; // mass of carbon,unit: kg
    scalar avogadro = 6.024e+23; // avogadro constant,unit: 1/mol
    scalar wC2H2 = 26.03824e-3;  // molar weight of C2H2, unit: kg/mol

    // read the concentrations for the species needed by the soot model (Mass fractions)
    const volScalarField& Y_C2H2 = mesh_.lookupObject<volScalarField>("C2H2");
    const volScalarField& Y_O2 = mesh_.lookupObject<volScalarField>("O2");
    const volScalarField& Y_OH = mesh_.lookupObject<volScalarField>("OH");
    const volScalarField& Y_CO = mesh_.lookupObject<volScalarField>("CO");
    const volScalarField& Y_H = mesh_.lookupObject<volScalarField>("H");
    const volScalarField& Y_H2 = mesh_.lookupObject<volScalarField>("H2");
    const volScalarField& Y_H2O = mesh_.lookupObject<volScalarField>("H2O");

    // determine the soot sub-model rates
    forAll(T,celli) // for loop in all the mesh(cell)
    {
    scalar Ti = T[celli];

    // real moments used for source term evaluation
    scalar realM0 = M0_[celli]*rho_[celli];
    scalar realM1 = M1_[celli]*rho_[celli];
    scalar realM2 = M2_[celli]*rho_[celli];
    scalar realM3 = M3_[celli]*rho_[celli];
    scalar realM4 = M4_[celli]*rho_[celli];

    // determine the inception rates
    if (Y_C2H2[celli]>0.0)
    {
    incepRateM0_[celli] = 2.2*sqrt(4.0*constant::mathematical::pi*kb/cMass/cperPAH_)
                          *pow(diamPAH_,2.0)*pow(avogadro,2)*sqrt(Ti)*
                          pow(Y_C2H2[celli]*rho_[celli]/wC2H2,2.0);
    }
    else
    {
      incepRateM0_[celli] = 0.0;
    }
    incepRateM1_[celli] = incepRateM0_[celli] * cperDimer_;
    incepRateM2_[celli] = incepRateM1_[celli] * cperDimer_;
    incepRateM3_[celli] = incepRateM2_[celli] * cperDimer_;
    incepRateM4_[celli] = incepRateM3_[celli] * cperDimer_;
    }

    // read the turbulence model to determine the effective eddy viscosity used in the soot model
    const compressible::turbulenceModel&
    turbModel = rho_.db().lookupObject<compressible::turbulenceModel>
    (
      turbulenceModel::propertiesName
    );

    // solve the soot moment transport equations
    // Transport equation for M0
    tmp<fvScalarMatrix> M0Eqn
    (
      fvm::ddt(rho_, M0_)
      +fvm::div(phi_, M0_)
      -fvm::laplacian(turbModel.muEff()/Sct_, M0_)
      ==
      incepRateM0_
    );
    solve(M0Eqn);

    // Transport equation for M1
    tmp<fvScalarMatrix> M1Eqn
    (
      fvm::ddt(rho_, M1_)
      +fvm::div(phi_, M1_)
      -fvm::laplacian(turbModel.muEff()/Sct_, M1_)
      ==
      incepRateM1_
    );
    solve(M1Eqn);

    // Transport equation for M2
    tmp<fvScalarMatrix> M2Eqn
    (
      fvm::ddt(rho_, M2_)
      +fvm::div(phi_, M2_)
      -fvm::laplacian(turbModel.muEff()/Sct_, M2_)
      ==
      incepRateM2_
    );
    solve(M2Eqn);

    // Transport equation for M3
    tmp<fvScalarMatrix> M3Eqn
    (
      fvm::ddt(rho_, M3_)
      +fvm::div(phi_, M3_)
      -fvm::laplacian(turbModel.muEff()/Sct_, M3_)
      ==
      incepRateM3_
    );
    solve(M3Eqn);

    // Transport equation for M4
    tmp<fvScalarMatrix> M4Eqn
    (
      fvm::ddt(rho_, M4_)
      +fvm::div(phi_, M4_)
      -fvm::laplacian(turbModel.muEff()/Sct_, M4_)
      ==
      incepRateM4_
    );
    solve(M4Eqn);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
