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
    along wi://mail.google.com/mail/u/0/#inboxth OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ThermalTopo1

Description
    Steady-state incompressible laminar flow solver with adjoint-based topology optimization
    for ducted flows. Modified from adjointshapeoptimizationFOAM to include temperature field and adjoint temperature field for solving multi-physics problem.
    optimization problem is to minimize the total pressure loss in the duct, while maximizing the heat transfer rate. 

    References:
  Fluid-thermal topology optimization of gas turbine blade internal cooling ducts
  S Ghosh, E Fernandez, J Kapat
  Journal of Mechanical Design 144 (5), 051703

  Topology Optimization of Internal Cooling Channels in Turbulent Flows
  S Ghosh
  University of Central Florida


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "wallFvPatch.H"
template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        laminarTransport.lookup("lambda") >> lambda;
        
		volScalarField alphat = turbulence->nut()/Prt;
       //alphat.correctBoundaryConditions();

            volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

      
         volScalarField gamma  = Foam::pow((1-(alpha/alphaMax)),2);
		 volScalarField Dac   = ((gamma*(alphaEff - DTc)) + DTc) ;	
		  volScalarField ass  = (fvc::grad(T)& U);
          // volScalarField head  = (fvc::grad(T)& fvc::grad(Ta));
            volScalarField head  = (alphaEff - DTc)*(fvc::laplacian(Ta,T));
        //alpha +=<F12>
        //    mesh.relaxationFactor("alpha")
        //   *(lambda*max(Ua & U, zeroSensitivity) - alpha);
        alpha +=
            mesh.fieldRelaxationFactor("alpha")
           *(min(max(alpha - lambda*(((Ua & U)-2/alphaMax*(1-(alpha/alphaMax))*(ass*Ta - head))+(alphaMax-2*alpha)*cons), zeroAlpha), alphaMax) - alpha);        //ass*Ta +

        zeroCells(alpha, inletCells);
        //zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)  
             ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

      



            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);



                // Temperature solver

     fvScalarMatrix TEqn
          (
                 fvm::ddt(T)
                + gamma*fvm::div(phi,T)
               == fvm::laplacian(Dac,T) 
           );

      TEqn.relax();  
     TEqn.solve();
	TEqn.relax();


            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
				
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }


           

       // Temperature solver

//   fvScalarMatrix TEqn
//        (
//               fvm::ddt(T)
//              + fvm::div(phi,T)
//             -fvm::laplacian(DTa,T) + alpha*(T-Twall)
//         );

//      TEqn.solve();
          
      
   //fvScalarMatrix Ta1Eqn
    //      (
     //          fvm::ddt(Ta)
      //           + fvm::div(-phi,Ta)
      //         -fvm::laplacian(Dac,Ta) 
      //     );
  // Ta1Eqn.relax();
   //    Ta1Eqn.solve();

           


     {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            //volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*fvc::dotInterpolate(fvc::snGrad(Ua), U)
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
              + fvm::Sp(alpha, Ua) + 1.2*1005*Ta*fvc::grad(T) 
             ==
                fvOptions(Ua)
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvOptions.constrain(UaEqn);

            solve(UaEqn == -fvc::grad(pa)) ;            

            fvOptions.correct(Ua);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);

      

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
            }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);
        }


       


        laminarTransport.correct();
        turbulence->correct();
        
          //adjoint Ta predictor

             fvScalarMatrix TaEqn
          (
                 fvm::ddt(Ta)
              + gamma* fvm::div(-phi,Ta)
              ==fvm::laplacian(Dac,Ta) 
           );
        TaEqn.relax();
        TaEqn.solve();
        sens = (Ua&U) ;  
        

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    
	
	Info<< "Starting Final calculation, Filtering intermediate porosities" << endl;
	
	forAll( mesh.C(), celli)
        {
				if (alpha[celli]>=20)
				{
				alpha[celli] = 200;
				}
				else 
				{
				alpha[celli] = 0;
				}
         		
		    
		}
		
		
	
	
	
	while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

       
		volScalarField alphat1 = turbulence->nut()/Prt;
        volScalarField gamma  = Foam::pow((1-(alpha/alphaMax)),2);
		
		
		 volScalarField alphaEff1("alphaEff", turbulence->nu()/Pr + alphat1);

		 volScalarField Dac   = ((gamma*(alphaEff1 - DTc)) + DTc) ;	

		
       //alphat.correctBoundaryConditions();

      
         

        zeroCells(alpha, inletCells);
        //zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)  
             ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

      



            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);



                // Temperature solver

     fvScalarMatrix TEqn
          (
                 fvm::ddt(T)
                + gamma*fvm::div(phi,T)
               == fvm::laplacian(Dac,T) 
           );

      TEqn.relax();  
      TEqn.solve();
	  TEqn.relax();


            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
				
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }


           

      


       


        laminarTransport.correct();
        turbulence->correct();
        
          

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }
	
	Info<< "End\n" << endl;

	
	
	
	
	
	
    return 0;
}


//************************************************************************ //
