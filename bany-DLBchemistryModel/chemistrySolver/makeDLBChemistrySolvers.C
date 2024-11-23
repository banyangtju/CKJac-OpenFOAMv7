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

#include "makeDLBChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeDLBChemistrySolverTypes(psiReactionThermo, constGasHThermoPhysics);
    makeDLBChemistrySolverTypes(psiReactionThermo, gasHThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeDLBChemistrySolverTypes(psiReactionThermo, icoPoly8HThermoPhysics);
    makeDLBChemistrySolverTypes(psiReactionThermo, constFluidHThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeDLBChemistrySolverTypes(psiReactionThermo, constHThermoPhysics);

    makeDLBChemistrySolverTypes(rhoReactionThermo, constGasHThermoPhysics);
    makeDLBChemistrySolverTypes(rhoReactionThermo, gasHThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeDLBChemistrySolverTypes(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeDLBChemistrySolverTypes(rhoReactionThermo, constFluidHThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeDLBChemistrySolverTypes(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeDLBChemistrySolverTypes(psiReactionThermo, constGasEThermoPhysics);
    makeDLBChemistrySolverTypes(psiReactionThermo, gasEThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeDLBChemistrySolverTypes(psiReactionThermo, icoPoly8EThermoPhysics);
    makeDLBChemistrySolverTypes(psiReactionThermo, constFluidEThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeDLBChemistrySolverTypes(psiReactionThermo, constEThermoPhysics);

    makeDLBChemistrySolverTypes(rhoReactionThermo, constGasEThermoPhysics);
    makeDLBChemistrySolverTypes(rhoReactionThermo, gasEThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeDLBChemistrySolverTypes(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeDLBChemistrySolverTypes(rhoReactionThermo, constFluidEThermoPhysics);
    makeDLBChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeDLBChemistrySolverTypes(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //