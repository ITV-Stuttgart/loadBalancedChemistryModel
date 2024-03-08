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
#include "pointToPointBuffer.H"
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
    public:
    // Struct to store the information required for sending cells to other
    // processors
    struct sendDataStruct
    {
        public:
            // processor to send information to
            label toProc=-1;

            // percent of own cells to send to processor toProc
            scalar percToSend=-1;
        
        // Constructor

            sendDataStruct() = default;

            sendDataStruct(scalar perc, label procId)
            : toProc(procId),percToSend(perc)
            {}
    };



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

        //- Store processor balancing data, first entry is a list of sendData
        //  information, second entry is a list of processor IDs of which 
        //  data is received
        Tuple2
        <
            List<sendDataStruct>,
            List<label>
        > 
        sendAndReceiveData_;

        //- List of size Pstream::nProcs() which is true if data is send to
        List<bool> sendToProcessor_;

        //- List of size Pstream::nProcs() which is true if data is received 
        //- from processor
        List<bool> receiveFromProcessor_;

        //- Needs to keep the buffer alive to not loose the information
        //  about send and receive sizes
        pointToPointBuffer pBufs_;

        //- Index of current iteration until the send/recv cell information
        //  is updated
        label iter_;

        //- Number of iterations till send/recv is updated
        label maxIterUpdate_;

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
        void updateProcessorBalancing();

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



