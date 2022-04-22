/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    calcTanalyticalSuckingProblem

Description
    Calculates and writes the analytical solution 
	for temperature field for one dimensional
	sucking problem. 
	It is assumed that the interface moves in the x direction.
Reference:
	@article{IRFAN2017132,
title = {A front tracking method for direct numerical simulation of evaporation process in a multiphase system},
journal = {Journal of Computational Physics},
volume = {337},
pages = {132-153},
year = {2017},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2017.02.036},
url = {https://www.sciencedirect.com/science/article/pii/S0021999117301304},
author = {Muhammad Irfan and Metin Muradoglu},
keywords = {Evaporation, Phase change, Front-tracking method, Multi-phase flows, The Clausius–Clapeyron relation, One-field formulation},
abstract = {A front-tracking method is developed for the direct numerical simulation of evaporation process in a liquid–gas multiphase system. One-field formulation is used to solve the flow, energy and species equations in the framework of the front tracking method, with suitable jump conditions at the interface. Both phases are assumed to be incompressible; however, the divergence-free velocity field condition is modified to account for the phase-change/mass-transfer at the interface. Both temperature and species gradient driven evaporation/phase-change processes are simulated. For the species gradient driven phase change process, the Clausius–Clapeyron equilibrium relation is used to find the vapor mass fraction and subsequently the evaporation mass flux at the interface. A number of benchmark cases are first studied to validate the implementation. The numerical results are found to be in excellent agreement with the analytical solutions for all the studied cases. The methods are then applied to study the evaporation of a static as well as a single and two droplets systems falling in the gravitational field. The methods are demonstrated to be grid convergent and the mass is globally conserved during the phase change process for both the static and moving droplet cases.}
}
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// modified from  wallHeatFlux
#include "singlePhaseTransportModel.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar func(scalar x, const scalar constValue)
{
	return x*Foam::exp(pow(x,2))*Foam::erf(x) - constValue;
}

scalar derivErf(const scalar x)
{
	return 2*Foam::exp(-pow(x,2))/Foam::sqrt(Foam::constant::mathematical::pi);
}

scalar derivFunc(scalar x)
{
	return 
		Foam::exp(pow(x,2))*Foam::erf(x) 
	  + x*(2*x*Foam::exp(pow(x,2))*Foam::erf(x) + Foam::exp(pow(x,2))*derivErf(x) ); 
}

// Interface position
dimensionedScalar X
       (
	       const scalar eps, 
		   const dimensionedScalar thermCond, 
		   const dimensionedScalar rho, 
		   const dimensionedScalar cp, 
		   const scalar time
	   )
{
	return 2*eps*Foam::sqrt(thermCond*time/rho/cp);
}

int main(int argc, char *argv[])
{

	argList::validArgs.append("time");
//    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    //instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"
#       include "createFields.H"
	#include "readTransportProperties.H"

	const scalar atime(readScalar(IStringStream(args.args()[1])()));

	const label maxIterNr = 20000;
	label iter = 0;
	const scalar tol = 1e-12;
	dimensionedScalar LHS, SquareBracketTerm;
	const dimensionedScalar alpha1(k1/rho1/cp1);
	const dimensionedScalar alpha2(k2/rho2/cp2);

	scalar epsilon = 0.1; // guess of starting value
	if (epsilonStartingValue)  
	{
	    Info << "Type the starting value for epsilon:" << endl;	
	    std::cin >> epsilon;
	}
	scalar epsilonPrev = 0;

	if (phaseChangeType == "evaporation")
	{
		SquareBracketTerm = ((Tinf-TSat)*cp2*k1*Foam::sqrt(alpha2)
			*Foam::exp(-epsilon*epsilon*rho2*rho2*alpha2/rho1/rho1/alpha1))
				/(hEvap*k2*Foam::sqrt(Foam::constant::mathematical::pi*alpha1)
					*Foam::erfc(epsilon*rho2*Foam::sqrt(alpha2)/rho1/Foam::sqrt(alpha1)));
		LHS = cp2*(Tw - TSat)/hEvap/Foam::sqrt(Foam::constant::mathematical::pi) 
			+ SquareBracketTerm*Foam::exp(epsilon*epsilon)*Foam::erf(epsilon);
	}
	// condensation not implemented
	//else if (phaseChangeType == "condensation")
	//{
	//	LHS = cp1*(TSat - Tw)/hEvap/Foam::sqrt(Foam::constant::mathematical::pi);
	//}
	else
	{
		FatalErrorIn("calcTanalyticalSuckingProblem") 
			//<< "phaseChangeType can be \"evaporation\" or \"condensation\" only"
			<< "phaseChangeType can be \"evaporation\" only (condensation not implemented so far)"
			<< nl
			<< exit(FatalError);
	}

	
	do 
	{
		epsilonPrev = epsilon;		
		SquareBracketTerm = ((Tinf-TSat)*cp2*k1*Foam::sqrt(alpha2)
			*Foam::exp(-epsilon*epsilon*rho2*rho2*alpha2/rho1/rho1/alpha1))
				/(hEvap*k2*Foam::sqrt(Foam::constant::mathematical::pi*alpha1)
					*Foam::erfc(epsilon*rho2*Foam::sqrt(alpha2)/rho1/Foam::sqrt(alpha1)));
		LHS = cp2*(Tw - TSat)/hEvap/Foam::sqrt(Foam::constant::mathematical::pi) 
			+ SquareBracketTerm*Foam::exp(epsilon*epsilon)*Foam::erf(epsilon);
		epsilon -= func(epsilon, LHS.value())/derivFunc(epsilon);	
		iter++;
	} while ( mag(epsilon - epsilonPrev) > tol && iter < maxIterNr );

    Info<< "Number of iterations: " << iter << endl;
	Info<< "epsilon = " << epsilon << endl;

	OFstream Tfile("T_SuckingProblem.txt");
	Tfile << "Time [s]\t"  << "x [m]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

	if (phaseChangeType == "evaporation")
	{
		volScalarField alphalInitial
		(
    	    IOobject
    	    (
				"alpha." + mixture->phase1Name(),
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::MUST_READ,
    	        IOobject::NO_WRITE
    	    ),
    	    mesh
		);

		dimensionedScalar avTau = k2/rho2/cp2*atime;
		dimensionedScalar alTau = k1/rho1/cp1*atime;
		dimensionedScalar analyticalInterfacePosition = X(epsilon,k2, rho2, cp2, atime);

	   	Info<< "\nAnalytical postion of the interface at time " 
			<< atime << " s is: " << analyticalInterfacePosition.value() << " m.\n"
			<< endl;
		
	   	Info<< "\nSaving the results for analytical temperature "
			<< "and alphal initial distribution for sucking problem\n" 
			<< "at time " << atime << " s" << endl;

		forAll(T, celli)
		{
			x[celli] = T.mesh().C()[celli].component(vector::X);
			if (x[celli] <= analyticalInterfacePosition.value())
			{
				T[celli] = TSat.value();
				alphalInitial[celli] = 0;
			}
			else
			{
				T[celli] = 
					Tinf.value() - (Tinf.value()-Tw.value())
					/(Foam::erfc(epsilon*rho2.value()*Foam::sqrt(alpha2.value())/rho1.value()/Foam::sqrt(alpha1.value()))) 
					*Foam::erfc(x[celli]/2.0/Foam::sqrt(alTau.value()) 
							+ epsilon*(rho2.value()-rho1.value())/rho1.value()*Foam::sqrt(alpha2.value()/alpha1.value()));

				alphalInitial[celli] = 1.0;
			}
		}
		T.write();
		alphalInitial.write();
	}

	//if (phaseChangeType == "condensation")
	//{
	//	volScalarField Tnum
	//	(
	//	    IOobject
	//	    (
	//	        "T",
	//	        runTime.timeName(),
	//	        mesh,
	//	        IOobject::MUST_READ,
	//	        IOobject::NO_WRITE
	//	    ),
	//	    mesh
	//	);

	//	dimensionedScalar tau = runTime.value();
	//	dimensionedScalar alTau = k1/rho1/cp1*tau;
	//	dimensionedScalar analyticalInterfacePosition= X(epsilon,k1, rho1, cp1, runTime.value());
	//	
	//   	Info<< "\nSaving the results for temperature to T.txt\n" << endl;
	//	forAll(analyticalTemperature, celli)
	//	{
	//		x[celli] = analyticalTemperature.mesh().C()[celli].component(vector::X);
	//		if (x[celli] > analyticalInterfacePosition.value())
	//		{
	//			analyticalTemperature[celli] = TSat.value();
	//		}
	//		else
	//		{
	//			analyticalTemperature[celli] = Tw.value() + (Tw.value() - TSat.value())
	//				/Foam::erf(epsilon)*Foam::erf(x[celli]/2.0/Foam::sqrt(alTau.value())); 
	//		}
	//
	//		Tfile << runTime.timeName() 
	//  		     << "\t" 
	//  			 << x[celli] 
	//  		     << "\t" 
	//  			 << Tnum[celli] 
	//  		     << "\t" 
	//  			 << analyticalTemperature[celli]
	//  		     << "\t" 
	//			 << mag(Tnum[celli] - analyticalTemperature[celli])/
	//			 		analyticalTemperature[celli]*100
	//  			 << endl;
	//	}
	//}

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
