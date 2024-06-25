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
    Extends the TDACChemistryModel with load balancing of the reation.

SourceFiles
    LoadBalancedTDACChemistryModel.H
    LoadBalancedTDACChemistryModel.C
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024

\*---------------------------------------------------------------------------*/
 
#ifndef LoadBalancedTDACChemistryModel_H
#define LoadBalancedTDACChemistryModel_H
 
#include "TDACChemistryModel.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "TDACDataContainer.H"
#include "pointToPointBuffer.H"
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
    public TDACChemistryModel<ReactionThermo, ThermoType>
{
public:

    // Struct to store the information required for sending cells to other
    // processors
    struct sendDataStruct
    {
        public:
            // processor to send information to
            label toProc=-1;

            // CPU time to send to processor
            scalar cpuTimeToSend=-1;
        
        // Constructor

            sendDataStruct() = default;

            sendDataStruct(scalar cpuTimeToSend, label procId)
            : toProc(procId),cpuTimeToSend(cpuTimeToSend)
            {}
    };


private:

    // Private Memeber Variables for TDAC

        //- Store for cpuStatistics
        clockTime clockTime_;
        scalar reduceMechCpuTime_;
        scalar addNewLeafCpuTime_;
        scalar growCpuTime_;
        scalar solveChemistryCpuTime_;
        scalar addToTableCpuTime_;
        scalar searchISATCpuTime_;
        scalar totalCpuTime_;

        //- Variable to store the solution vector. 
        //  Allocated here to avoid reallocation
        scalarField Rphiq_;

    // Private Member Variables for Load Balancing

        //- List of TDACDataContainer for each cell of the field
        //  Always has the size of the field
        DynamicList<TDACDataContainer> cellDataField_;

        //- Switch to check if it is called the first time in the simulation
        bool firstTime_{true};

        //- Number of cells that were not found in the table and need to be
        //  solved
        label cellsToSolve_;

        //- Store the number of species of the reduced mechanism
        label nSpecieReduced_;

        //- Store processor balancing data, first entry is a list of sendData
        //  information, second entry is the recv processor ID
        Tuple2
        <
            List<sendDataStruct>,
            List<label>
        > 
        sendAndReceiveData_;

        //- Store remote computed cells in a separate table
        autoPtr<chemistryTabulationMethod<ReactionThermo, ThermoType>>
            tabulationRemote_;

        //- List of size Pstream::nProcs() which is true if data is send to
        List<bool> sendToProcessor_;

        //- List of size Pstream::nProcs() which is true if data is received 
        //- from processor
        List<bool> receiveFromProcessor_;

        //- Needs to keep the buffer alive to not loose the information
        //  about send and receive sizes
        pointToPointBuffer pBufs_;

        //- Index of current iteration before the proc-tp-proc connection
        //  is recalculated for the load balancing
        label iter_;

        //- Maximum number of iterations before the proc-to-proc load balancing
        //  is updated
        label maxIterUpdate_;

    // Private Member Functions

        //- Add cell to cellDataList for parallel processing
        void addCell
        (
            DynamicList<TDACDataContainer*>& cellList,
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
            const DynamicList<TDACDataContainer*>& cellList,
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
        //  List[index].second: Returns a List<scalar> with:
        //                        List[0]: processor ID
        //                        List[1]: Number of cells to compute 
        //                        List[2]: Time to search ISAT table
        List<std::pair<scalar,List<scalar>>> 
        getSortedCPUTimesOnProcessor() const;
        
        //- Update the totalCpuTime_ variable 
        void updateTotalCpuTime
        (
            const DynamicList<TDACDataContainer*>& cellList
        );

        //- Add cell to ISAT table -- after solving
        //  Switch to set if local or remote cells are solved
        //  Switch if the reduced mechanism needs to be recomputed
        //  this is required if the cell was sent between processors
        void addCellToTable
        (
            const TDACDataContainer& cData,
            const bool isLocal,
            const bool requiresRecomputeReducedMech=true
        );

        //- Lookup the cell data in the ISAT table
        bool lookUpCellInTable
        (
            TDACDataContainer& cData,
            const bool isLocal
        );

        //- Solve the reaction for all cells in the given list
        //  Flag sets if it is local or remote cell computation
        void solveCellList
        (
            UList<TDACDataContainer*>& cellList,
            const bool isLocal
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
        LoadBalancedTDACChemistryModel
        (
            const LoadBalancedTDACChemistryModel&
        ) = delete;
  
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

    // ODE functions (overriding functions in TDACChemistryModel to take

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



