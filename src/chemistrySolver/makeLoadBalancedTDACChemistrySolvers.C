/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "makeLoadBalancedTDACChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constGasHThermoPhysics);
    
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, gasHThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    )
    ;
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, icoPoly8HThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constFluidHThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constHThermoPhysics);


    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constGasHThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, gasHThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constFluidHThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constHThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constGasEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, gasEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, icoPoly8EThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constFluidEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(psiReactionThermo, constEThermoPhysics);

    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constGasEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, gasEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, icoPoly8EThermoPhysics);

    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constFluidEThermoPhysics);
    makeLoadBalancedTDACChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeLoadBalancedTDACChemistrySolverTypes(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
