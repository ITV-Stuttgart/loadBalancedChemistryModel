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

Class
    Foam::TDACDataContainer
 
Description
    Data container to store the species, pressure, and temperature information
    of cells for processor streaming with Pstream.

SourceFiles
    TDACDataContainer.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef TDACDataContainer_H
#define TDACDataContainer_H

#include "fvCFD.H"
#include "Ostream.H"
#include "Istream.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                IO Functions required for Parallel Handling
\*---------------------------------------------------------------------------*/

class TDACDataContainer;

// Read and write the eulerian fields data class
Istream& operator >> (Istream& is, TDACDataContainer& eField);

Ostream& operator << (Ostream& os, const TDACDataContainer& eField);

// When comparing TDACDataContainer always return untrue as otherwise all 
// fields have to be checked
// --> This operation is only needed to write out a list of TDACDataContainer
//     which checks for uniform list entries
inline bool operator==(const TDACDataContainer& e1, const TDACDataContainer& e2)
{
    return false;
}

inline bool operator!=(const TDACDataContainer& e1, const TDACDataContainer& e2)
{
    return true;
}


/*---------------------------------------------------------------------------*\
                     Class TDACDataContainer
\*---------------------------------------------------------------------------*/
// Base struct to store the eulerian field properties for each particle

class TDACDataContainer
{
    protected:

        //- Composition vector (Y, T, p)
        Field<scalar> phiq_;

        //- Concentration
        Field<scalar> c_;

        //- Concentration prior solving
        Field<scalar> c0_;

        //- Temperature
        scalar T_;

        //- Pressure
        scalar p_;

        //- Density
        scalar rho_;

        //- Time deltaT
        scalar deltaT_;

        //- Chemical time scale 
        scalar deltaTChem_;

        //- cpu time required for computation - default init to zero
        scalar cpuTime_{0};

        //- Store processor ID
        label procID_;

        //- Label of cellID
        label cellID_;

        //- Is it a local particle
        //  set to false if transferred
        bool local_{true};
    
    public:
    
        TDACDataContainer() = default;

        //- Construct with given size for lists 
        TDACDataContainer
        (
            const label nSpecies,
            const label nAdditions
        )
        :
            phiq_(nSpecies+nAdditions),
            c_(nSpecies),
            c0_(nSpecies)
        {};
        
        //- Construct from input stream
        TDACDataContainer(Istream& is);
        
        // Functions for IO operations
        
            //- Read TDACDataContainer from stream
            Istream& read(Istream& is);
            
            //- Write TDACDataContainer to stream
            Ostream& write(Ostream& os) const;
           

       // Modify

            //- Return species
            Field<scalar>& phiq() {return phiq_;}

            //- Return concentration
            Field<scalar>& c() {return c_;}

            //- Return old concentration
            Field<scalar>& c0() {return c0_;}

            //- Access the temperature information
            scalar& T() {return T_;}

            //- Access the pressure information
            scalar& p() {return p_;}

            //- Access the density
            scalar& rho() {return rho_;}

            //- Access deltaT
            scalar& deltaT() {return deltaT_;}

            //- Access deltaTChem() 
            scalar& deltaTChem() {return deltaTChem_;}

            //- Access processor ID
            label& proc() {return procID_;}

            //- Access cpu time
            scalar& cpuTime() {return cpuTime_;}

            //- Acess cellID
            label& cellID() {return cellID_;}

    // Const Access

            //- Return species
            const Field<scalar>& phiq() const {return phiq_;}

            //- Return concentration
            const Field<scalar>& c() const {return c_;}

            //- Return old concentration
            const Field<scalar>& c0() const {return c0_;}

            //- Access the temperature information
            const scalar& T() const {return T_;}

            //- Access the pressure information
            const scalar& p() const {return p_;}

            //- Access the density
            const scalar& rho() const {return rho_;}

            //- Access deltaT
            const scalar& deltaT() const {return deltaT_;}

            //- Access deltaTChem() 
            const scalar& deltaTChem() const {return deltaTChem_;}

            //- Access processor ID
            const label& proc() const {return procID_;}

            //- Access cpu time
            const scalar& cpuTime() const {return cpuTime_;}

            //- Acess cellID
            const label& cellID() const {return cellID_;}

    // Check if container had been streamed 

            //- check if the container is local
            inline bool& local() {return local_;}
        
            //- Const access if it is a local container
            inline const bool& local() const {return local_;}
};

}   // End of namespace Foam
#endif
