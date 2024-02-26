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
    Foam::LoadBalancedChemistryModel

Description
    Extends the StandardChemistryModel with load balancing of the reation.

SourceFiles
    LoadBalancedChemistryModel.H
    LoadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/
 
#ifndef LoadBalancedChemistryModel_H
#define LoadBalancedChemistryModel_H
 
#include "StandardChemistryModel.H"
#include "baseDataContainer.H"
#include "OFstream.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                      Class LoadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/
 
template<class ReactionThermo, class ThermoType>
class LoadBalancedChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
private:

    // Private Member Variables

        //- Minimum percent of cells to send to processors
        //  Only processors that would receive more the n-percent of cells 
        //  are considered for load balancing
        //  --> Value is given as absolute fraction, not in percent!
        //      Default value is: 0.02
        scalar minFractionOfCellsToSend_;

        //- List to store the cell information in the baseDataContainer
        DynamicList<baseDataContainer> cellDataList_;

        //- Get list of cells on each processor
        List<label> cellsOnProcessors_;

        //- Switch to check if it is called the first time in the simulation
        bool firstTime_{true};

        //- Total cpu time on all processors
        scalar totalCpuTime_;

    // Private Member Functions

        //- Build the baseDataContainer list from the cells
        template<class DeltaTType>
        void buildCellDataList(const DeltaTType&);
        
        //- Update teh cellDataList with new cell values
        template<class DeltaTType>
        void updateCellDataList(const DeltaTType&);

        //- Calculate the cells to send/recv to/from other 
        //- processors
        //  Returns a list of sublists of cells to send to processor i
        void cellsToSend
        (
            const DynamicList<baseDataContainer>& cellList,
            const scalar cpuTimeToSend,
            const label& start,
            label& end
        );
        
        //- Get percentage of processor data to send to other processors
        //- to balance the processor load
        Tuple2
        <
            List<Tuple2<scalar,label>>,
            List<label>
        > 
        getProcessorBalancing();

        //- Get cpu times of each processor in a descending order
        //  Format:
        //  List[index].first:  CPU time required for last time step
        //  List[index].second: Returns an Foam::Pair<label> with
        //                        Pair.first():  processor ID
        //                        Pair.second(): Number of particles on that processor 
        List<std::pair<scalar,Pair<label>>> getSortedCPUTimesOnProcessor() const;
        
        //- Update the totalCpuTime_ variable 
        void updateTotalCpuTime
        (
            const DynamicList<baseDataContainer>& reactParList
        );

        //- solve the reaction for all cells in the given list
        void solveCellList
        (
            UList<baseDataContainer>& cellList
        );

        //- Solve chemistry for once cell
        void solveCell(baseDataContainer& cellData);


        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

public:
    //- Runtime type information
    TypeName("LoadBalancedChemistryModel");

    // Constructors

        //- No copy construct
        LoadBalancedChemistryModel(const LoadBalancedChemistryModel&) = delete;
  
        //- No copy assignment
        void operator=(const LoadBalancedChemistryModel&) = delete;
  
        //- Construct from thermo
        LoadBalancedChemistryModel(ReactionThermo& thermo);
  
  
        //- Destructor
        virtual ~LoadBalancedChemistryModel() {};
        
        //- Solve the reaction system for the given time step with load
        //  balancing and return the characteristic time
        virtual scalar solve(const scalar deltaT);

        //- Solve the reaction system for the given time step with load
        //  balancing and return the characteristic time
        virtual scalar solve(const scalarField& deltaT);


        //- Access the stored cell data
        DynamicList<baseDataContainer>& cellDataList()
        {
            return cellDataList_;
        }

        //- Const access the stored cell data
        const DynamicList<baseDataContainer>& cellDataList() const
        {
            return cellDataList_;
        }

    // ODE functions (overriding abstract functions in ODE.H)
 
        virtual void solve
        (
            scalarField &c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const = 0;
};

}   // End of namespace Foam

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
 #ifdef NoRepository
     #include "LoadBalancedChemistryModel.C"
 #endif

#endif



