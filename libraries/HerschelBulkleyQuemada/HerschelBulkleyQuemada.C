/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "HerschelBulkleyQuemada.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HerschelBulkleyQuemada, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HerschelBulkleyQuemada,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::HerschelBulkleyQuemada::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
    volScalarField epssStar = U_.mesh().lookupObject<volScalarField>("epss"); 

epssStar.min((epssMax_.value()-0.01)); 
epssStar.max(0.01); 

volScalarField Tau0 = tau0_*pow( (epssStar/epssMax_) / max((1-epssStar/epssMax_),SMALL ),  Dstar_ ); 

volScalarField nufInf = (nuf_/pow(max(SMALL,(1-epssStar/epssMax_)),2));

return
    (
        min
        (
            nu0_,
	    (Tau0 + nufInf*rtone*pow(tone*sr(), n_))
           /(max(sr(), dimensionedScalar ("small", dimless/dimTime, SMALL)))
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::HerschelBulkleyQuemada::HerschelBulkleyQuemada
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HerschelBulkleyQuemadaCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    nuf_("nuf", dimViscosity, HerschelBulkleyQuemadaCoeffs_),
    epssMax_("epssMax", dimless, HerschelBulkleyQuemadaCoeffs_),
    n_("n", dimless, HerschelBulkleyQuemadaCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, HerschelBulkleyQuemadaCoeffs_),
    Dstar_("Dstar", dimless, HerschelBulkleyQuemadaCoeffs_),
    nu0_("nu0", dimViscosity, HerschelBulkleyQuemadaCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::HerschelBulkleyQuemada::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HerschelBulkleyQuemadaCoeffs_ =
        viscosityProperties.subDict(typeName + "Coeffs");

    HerschelBulkleyQuemadaCoeffs_.lookup("nuf") >> nuf_;
    HerschelBulkleyQuemadaCoeffs_.lookup("epssMax") >> epssMax_;
    HerschelBulkleyQuemadaCoeffs_.lookup("n") >> n_;
    HerschelBulkleyQuemadaCoeffs_.lookup("tau0") >> tau0_;
    HerschelBulkleyQuemadaCoeffs_.lookup("Dstar") >> Dstar_;
    HerschelBulkleyQuemadaCoeffs_.lookup("nu0") >> nu0_;

    return true;
}


// ************************************************************************* //

