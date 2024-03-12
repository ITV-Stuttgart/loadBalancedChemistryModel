/*---------------------------------------------------------------------------*\
  ==========   ==========   |
      ||           ||       | An OpenFOAM Extension
      ||  Tabular  ||       |
      ||  Thermo   ||       | Copyright (C) 2018 Jan Gärtner
      ||           ||       |
-------------------------------------------------------------------------------

License
    This file is part of tabularThermo.

    tabularThermo is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    tabularThermo is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with  tabularThermo.  If not, see <http://www.gnu.org/licenses/>.

Description
    Test the load-balanced standard chemistry model by comparison to the 
    standard model. 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/



// Includes for Catch2 unit testing framework
#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// For global arguments
#include "globalFoamArgs.H"

// Standard C++ includes
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <iostream>

// OpenFOAM includes
#include "fvCFD.H"
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "ode.H"
#include "LoadBalancedChemistryModel.H"


TEST_CASE("standardChemistryModel-Test","[chemistry]")
{
    // =========================================================================
    //                      Prepare Case
    // =========================================================================
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    #include "createMesh.H"

    // Create a thermo model
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    // Create a typedef
    using chemModelLB = 
        LoadBalancedChemistryModel<psiReactionThermo,gasHThermoPhysics>;

    using chemModelStd = 
        StandardChemistryModel<psiReactionThermo,gasHThermoPhysics>;

    // Create the load-balanced and standard chemistry model
    ode<chemModelLB> cModelLB(thermo);
    ode<chemModelStd> cModelStd(thermo);

    const scalar deltaT = 1E-6;

    cModelLB.chemModelLB::solve(deltaT);
    cModelStd.chemModelStd::solve(deltaT);

    // Compare the computed reaction rate
    for (label specieI=0; specieI < cModelLB.nSpecie(); specieI++)
    {
        auto RRLB = cModelLB.RR(specieI);
        auto RRStd = cModelStd.RR(specieI);
        forAll(RRLB,celli)
        {
            REQUIRE_THAT
            (
                RRLB[celli],
                Catch::Matchers::WithinRel(RRStd[celli],1E-6)
            );
        }
    }
}

