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
    the Free SoftwaRe Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SpearmanQuemada.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(SpearmanQuemada, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        SpearmanQuemada,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SpearmanQuemada::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
    volScalarField epssStar = U_.mesh().lookupObject<volScalarField>("epss"); 

//making sure that epssEffStar has maximim and minimum values for stability//

epssStar.min((epssMax_.value()-0.01)); //0.05
epssStar.max(0.3); //0.01

//Creating Tau and nuf fields
volScalarField Tau0 = tau0_*pow( (epssStar/epssMax_) / max((1-epssStar/epssMax_),SMALL ),  Dstar_ ); 
volScalarField nufInf = (nuf_/pow(max(SMALL,(1-epssStar/epssMax_)),2));  

//calculating the viscosity values
return
(
min(nu0_,
  pow(	
	(   
	    //First Term
	    pow(Tau0,(1/Dstar_))

	    //Second Term
	  + pow(
		1+
                ( 
		1/(
		   Beta_*max(
			    pow(sr(),m_) , pow(dimensionedScalar ("small", dimless/dimTime, SMALL),m_)
			    ) 
                  )
		)*
		pow(Tau0/nufInf,m_),
                (1/Dstar_))

	   //Third Term
	  * pow(nufInf*sr(),(1/Dstar_))
	)
    ,Dstar_)

       //Dividing by strain to get an effective viscosity
      /(max(sr(), dimensionedScalar ("small", dimless/dimTime, SMALL)))
) 
);
}

//where tau0_ is the critical stress constant and nuf_ is the viscosity of water
//Following Spearman paper:
//epssEffMax is the maximum effective solid fraction, D is a fractal dimension constant
//Dstar is 4-D


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::SpearmanQuemada::SpearmanQuemada
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    SpearmanQuemadaCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    nuf_("nuf", dimViscosity, SpearmanQuemadaCoeffs_),
    epssMax_("epssMax", dimless, SpearmanQuemadaCoeffs_),
    n_("n", dimless, SpearmanQuemadaCoeffs_),
    m_("m", dimless, SpearmanQuemadaCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, SpearmanQuemadaCoeffs_),
    Dstar_("Dstar", dimless, SpearmanQuemadaCoeffs_),
    nu0_("nu0", dimViscosity, SpearmanQuemadaCoeffs_),
    Beta_("Beta", dimless, SpearmanQuemadaCoeffs_),
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

bool Foam::viscosityModels::SpearmanQuemada::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    SpearmanQuemadaCoeffs_ =
        viscosityProperties.subDict(typeName + "Coeffs");

    SpearmanQuemadaCoeffs_.lookup("nuf") >> nuf_;
    SpearmanQuemadaCoeffs_.lookup("epssMax") >> epssMax_;
    SpearmanQuemadaCoeffs_.lookup("n") >> n_;
    SpearmanQuemadaCoeffs_.lookup("m") >> m_;
    SpearmanQuemadaCoeffs_.lookup("tau0") >> tau0_;
    SpearmanQuemadaCoeffs_.lookup("Dstar") >> Dstar_;
    SpearmanQuemadaCoeffs_.lookup("nu0") >> nu0_;
    SpearmanQuemadaCoeffs_.lookup("Beta") >> Beta_;

    return true;
}


// ************************************************************************* //

//other stuff
//min
        //(
           // nu0_,
            //(
            //max(
	//	tauMin_, tau0_*pow((epssEff/epssEffMax_)/(1-epssEff/epssEffMax_),Dstar_)
	 //      )
         //   + (k_/pow((1-epssEff/epssEffMax_),2))*rtone*pow(tone*sr(), n_)
         //   )
         //  /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
