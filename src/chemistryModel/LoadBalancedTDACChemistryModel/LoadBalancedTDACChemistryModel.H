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
    Foam::LoadBalancedTDACChemistryModel

Description
    Extends the StandardChemistryModel with load balancing of the reation.

SourceFiles
    LoadBalancedTDACChemistryModel.H
    LoadBalancedTDACChemistryModel.C

\*---------------------------------------------------------------------------*/
 
#ifndef LoadBalancedTDACChemistryModel_H
#define LoadBalancedTDACChemistryModel_H
 
#include "TDACChemistryModel.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "TDACDataContainer.H"
#include "OFstream.H"
#include "clockTime.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                      Class LoadBalancedTDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/
 
template<class ReactionThermo, class ThermoType>
class LoadBalancedTDACChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
private:

    // Private Memeber Variables for TDAC

        //- Store for cpuStatistics
        clockTime clockTime_;
        scalar reduceMechCpuTime_;
        scalar addNewLeafCpuTime_;
        scalar growCpuTime_;
        scalar solveChemistryCpuTime_;
        scalar searchISATCpuTime_;

    // Private Member Variables for Load Balancing

        //- List to store the cell information in the TDACDataContainer
        DynamicList<TDACDataContainer> cellDataList_;

        //- Switch to check if it is called the first time in the simulation
        bool firstTime_{true};

        //- Total cpu time on all processors
        scalar totalCpuTime_;

    // Private Member Functions

        //- Add cell to cellDataList for parallel processing
        void addCell
        (
            const scalarField& phiq,
            const scalar& T,
            const scalar& p,
            const scalar& rho,
            const scalar& deltaT,
            const label& celli
        );

        //- Calculate the cells to send/recv to/from other 
        //- processors
        //  Returns a list of sublists of cells to send to processor i
        void cellsToSend
        (
            const DynamicList<TDACDataContainer>& cellList,
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
            const DynamicList<TDACDataContainer>& reactParList
        );

        //- solve the reaction for all cells in the given list
        void solveCellList
        (
            UList<TDACDataContainer>& cellList
        );

        //- Solve chemistry for once cell
        void solveCell(TDACDataContainer& cellData);


        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

public:
    //- Runtime type information
    TypeName("LoadBalancedTDAC");

    // Constructors

        //- No copy construct
        LoadBalancedTDACChemistryModel(const LoadBalancedTDACChemistryModel&) = delete;
  
        //- No copy assignment
        void operator=(const LoadBalancedTDACChemistryModel&) = delete;
  
        //- Construct from thermo
        LoadBalancedTDACChemistryModel(ReactionThermo& thermo);
  
  
        //- Destructor
        virtual ~LoadBalancedTDACChemistryModel() {};
        
        //- Solve the reaction system for the given time step with load
        //  balancing and return the characteristic time
        virtual scalar solve(const scalar deltaT);

        //- Solve the reaction system for the given time step with load
        //  balancing and return the characteristic time
        virtual scalar solve(const scalarField& deltaT);


        //- Access the stored cell data
        DynamicList<TDACDataContainer>& cellDataList()
        {
            return cellDataList_;
        }

        //- Const access the stored cell data
        const DynamicList<TDACDataContainer>& cellDataList() const
        {
            return cellDataList_;
        }

    // ODE functions (overriding functions in StandardChemistryModel to take

        virtual void solve
        (
            scalarField& c,
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
    #include "LoadBalancedTDACChemistryModel.C"
#endif

#endif



