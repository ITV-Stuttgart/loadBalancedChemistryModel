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

#include "makeLoadBalancedChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    )
    ;
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constHThermoPhysics
    );


    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constHThermoPhysics
    );


    // Chemistry solvers based on sensibleInternalEnergy
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        psiReactionThermo,
        constEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeLoadBalancedChemistrySolverTypes
    (
        rhoReactionThermo,
        constEThermoPhysics
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
