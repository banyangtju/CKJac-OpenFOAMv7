/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "LoadBalancedChemistryModel.H"
#include "LoadBalanced_pyJacChemistryModel.H"
#include "LoadBalanced_CKJacChemistryModel.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


    // Make base types
//    makeChemistryModel(psiReactionThermo);
//    makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    
    
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );

    
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );

    
    


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );



    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalancedChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );



    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    ); 
    
    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

 
    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    
    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );

    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );



    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_pyJacChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    ); 
	
	
	
	
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    ); 
    
    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

 
    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    
    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );

    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );



    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        LoadBalanced_CKJacChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );	
}

// ************************************************************************* //