/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021,2023 OpenCFD Ltd.
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

#include "LoadBalancedChemistryModel.H"

template<class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::LoadBalancedChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
{
    cellsOnProcessors_.resize(Pstream::nProcs());
    // Gather the number of particles on each processor
    List<label> cellsOnProcessors(Pstream::nProcs());
    cellsOnProcessors[Pstream::myProcNo()] = this->mesh().C().size();
    
    Pstream::gatherList(cellsOnProcessors);
    Pstream::scatterList(cellsOnProcessors);

    cellsOnProcessors_ = cellsOnProcessors;
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
void Foam::LoadBalancedChemistryModel
<ReactionThermo, ThermoType>::buildCellDataList
(
    const DeltaTType& deltaT
)
{
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    // Reserve at least enough space for all cells on the local mesh
    cellDataList_.reserve(rho.size());

    forAll(rho, celli)
    {
        // store the cell data in a container
        baseDataContainer cData(this->nSpecie_);
        cellDataList_.append(cData);
    }

    updateCellDataList(deltaT);
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
void Foam::LoadBalancedChemistryModel
<ReactionThermo, ThermoType>::updateCellDataList
(
    const DeltaTType& deltaT
)
{
    label MyProcNo = Pstream::myProcNo();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    // Additional checks if FULLDEBUG is activated:
    #ifdef FULLDEBUG
        // Check that the number of cells has not changed
        if (rho.size() != cellDataList().size())
            FatalError 
                << "updateCellDataList does not work if the mesh changes"
                << "  -- mesh size: " << p.size() << " cellDataList size: "
                << cellDataList().size()
                << exit(FatalError);
    #endif

    scalarField c0(this->nSpecie_);

    // Loop over all entries
    forAll(cellDataList_,celli)
    {
        auto& cellData = cellDataList_[celli];

        cellData.proc() = MyProcNo;

        cellData.T() = T[celli];

        cellData.p() = p[celli];

        cellData.rho() = rho[celli];

        cellData.deltaT() = deltaT[celli];

        cellData.deltaTChem() = this->deltaTChem_[celli];

        cellData.cellID() = celli;

        // Set species
        forAll(this->Y_,j)
        {
            cellData.Y()[j] = this->Y_[j][celli];
        }
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::cellsToSend
(
    const DynamicList<baseDataContainer>& cellList,
    const scalar cpuTimeToSend,
    const label& start,
    label& end
)
{
    end = cellList.size();
    scalar cpuTime = 0;
    
    // go from start index and add as many particles until the cpuTimeToSend
    // is reached
    for (label i=start; i < cellList.size(); i++)
    {
        if (cpuTime >= cpuTimeToSend)
        {
            end = i;
            break;
        }

        cpuTime += cellList[i].cpuTime();
    }
}


template<class ReactionThermo, class ThermoType>
Foam::Tuple2
<
    Foam::List<Foam::Tuple2<scalar,label>>,
    Foam::List<label>
>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getProcessorBalancing()
{
    // Number of processors 
    const scalar numProcs = Pstream::nProcs();

    auto sortedCpuTimeOnProcessors = getSortedCPUTimesOnProcessor();

    // calculate average time spent on each cpu
    scalar averageCpuTime = 0;
    forAll(sortedCpuTimeOnProcessors,i)
    {
        averageCpuTime += sortedCpuTimeOnProcessors[i].first;
    }
    averageCpuTime /= numProcs;
    
    Info << "Average CPU time : "<<averageCpuTime<<endl;

    // list of the distributed load for all processors
    List<DynamicList<Tuple2<scalar,label>>> distributedLoadAllProcs(numProcs);
    
    // list of processors to receive data from
    List<DynamicList<label>> receiveDataFromProc(numProcs);

    // balance the load by calculating the percentages to be send 
    forAll(sortedCpuTimeOnProcessors,i)
    {
        const label procI = sortedCpuTimeOnProcessors[i].second.first();
        const label nCellsOnProcI = sortedCpuTimeOnProcessors[i].second.second();

        // List of processors to send information to
        DynamicList<Tuple2<scalar,label>>& sendLoadList = 
            distributedLoadAllProcs[procI];
        
        // Reserve space
        sendLoadList.reserve
        (
            std::floor(0.5*(numProcs-i))
        );
        
        
        const scalar cpuTimeProcI = sortedCpuTimeOnProcessors[i].first;
        
        scalar cpuTimeOverhead = cpuTimeProcI - averageCpuTime;
        
        
        // Loop over the other processors and distribute the load
        // in reverse order
        for (label k=numProcs-1; k > 0; k--)
        {
            const label procK = sortedCpuTimeOnProcessors[k].second.first();
            
            const scalar cpuTimeProcK = sortedCpuTimeOnProcessors[k].first;
            
            const scalar capacityOfProcK = averageCpuTime - cpuTimeProcK;
            
            if (capacityOfProcK <= 0)
                continue;
            
            const scalar newCapacity = capacityOfProcK - cpuTimeOverhead;
            
            // Check that the number of particles to send is greater than 1
            if ((cpuTimeOverhead/cpuTimeProcI*nCellsOnProcI) < 2)
                continue;
            

            if (newCapacity > 0)
            {                
                // add the cpuTimeOverhead to the load of this processor
                sortedCpuTimeOnProcessors[k].first += cpuTimeOverhead;
                sortedCpuTimeOnProcessors[i].first -= cpuTimeOverhead;
                
                // Processor procI sends information to procK
                // Therefore the receive data list of procK has to be updated
                receiveDataFromProc[procK].reserve(0.5*numProcs);
                receiveDataFromProc[procK].append(procI);

                sendLoadList.append
                (
                    Tuple2<scalar,label>
                    (
                        cpuTimeOverhead/cpuTimeProcI,
                        procK
                    )
                );

                cpuTimeOverhead = 0;

                break;
            }
            else
            {
                sortedCpuTimeOnProcessors[k].first += capacityOfProcK;
                sortedCpuTimeOnProcessors[i].first -= capacityOfProcK;
                cpuTimeOverhead -= capacityOfProcK;

                // Processor procI sends information to procK
                // Therefore the receive data list of procK has to be updated
                receiveDataFromProc[procK].reserve(0.5*numProcs);
                receiveDataFromProc[procK].append(procI);
                
                sendLoadList.append
                (
                    Tuple2<scalar,label>
                    (
                        capacityOfProcK/cpuTimeProcI,
                        sortedCpuTimeOnProcessors[k].second.first()
                    )
                );
            }
        }
    }
    
    // return only the list for the current processor
    return Tuple2
    <
        List<Foam::Tuple2<scalar,label>>,
        List<label>
    >
    (
        distributedLoadAllProcs[Pstream::myProcNo()],
        receiveDataFromProc[Pstream::myProcNo()]
    );
}


template<class ReactionThermo, class ThermoType>
Foam::List<std::pair<scalar,Foam::Pair<label>>> 
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getSortedCPUTimesOnProcessor() const
{
    // Number of processors 
    const scalar numProcs = Pstream::nProcs();
        
    // Gather the data from all processors
    List<scalar> cpuTimeOnProcessors(numProcs);
    cpuTimeOnProcessors[Pstream::myProcNo()] = totalCpuTime_;


    Pstream::gatherList(cpuTimeOnProcessors);
    Pstream::scatterList(cpuTimeOnProcessors);

    // use std::pair for std::sort algorithm
    List<std::pair<scalar,Pair<label>>> sortedCpuTimeOnProcessors(numProcs);
    
    forAll(sortedCpuTimeOnProcessors,i)
    {
        sortedCpuTimeOnProcessors[i].first = cpuTimeOnProcessors[i];
        sortedCpuTimeOnProcessors[i].second.first() = i;
        sortedCpuTimeOnProcessors[i].second.second() = cellsOnProcessors_[i];
    }
    
    // sort using std::sort()
    // first entry has largest cpu time --> descending order
    std::stable_sort
    (
        sortedCpuTimeOnProcessors.begin(), 
        sortedCpuTimeOnProcessors.end(),
        std::greater<std::pair<scalar,Pair<label>>>()
    );

    return sortedCpuTimeOnProcessors;
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateTotalCpuTime
(
    const DynamicList<baseDataContainer>& reactCellList
)
{
    // Calculate the total time spent solving the particles on this processor
    totalCpuTime_  = 0;
    
    for (const auto& obj : reactCellList)
        totalCpuTime_ += obj.cpuTime();
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveCell
(
    baseDataContainer& cellData
)
{
    if (cellData.T() > this->Treact_)
    {
        // We can use here the specieThermo at any processor
        // as the weight of a specie is static and the same on each processor
        auto& Y = cellData.Y();
        
        scalarField c0(this->nSpecie_);

        for (label i=0; i<this->nSpecie_; i++)
        {
            this->c_[i] = cellData.rho()*Y[i]/this->specieThermo_[i].W();
            c0[i] = this->c_[i];
        }

        // Initialise time progress
        scalar timeLeft = cellData.deltaT();

        // Calculate the chemical source terms
        while (timeLeft > SMALL)
        {
            scalar dt = timeLeft;
            this->solve(this->c_, cellData.T(), cellData.p(), dt, cellData.deltaTChem());
            timeLeft -= dt;
        }

        cellData.deltaTChem() =
            min(cellData.deltaTChem(), this->deltaTChemMax_);

        for (label i=0; i<this->nSpecie_; i++)
        {
            cellData.RR()[i] =
                (this->c_[i] - c0[i])*this->specieThermo_[i].W()/cellData.deltaT();
        }
    }
    else
    {
        for (label i=0; i<this->nSpecie_; i++)
        {
            cellData.RR()[i] = 0;
        }
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveCellList
(
    UList<baseDataContainer>& cellList
)
{    
    for(baseDataContainer& cellData : cellList)
    {
        // We cannot use here cpuTimeIncrement() of OpenFOAM as this 
        // returns only measurements in 100Hz or 1000Hz intervals depending
        // on the installed kernel 
        auto start = std::chrono::high_resolution_clock::now();

        solveCell
        (
            cellData
        );
        
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end-start));

        // Add the time required to solve this particle to the list 
        // as seconds
        cellData.cpuTime() = duration.count()*1.0E-6;
    }
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    // if it is not a prallel run exit with an error
    if (!Pstream::parRun())
        FatalError 
            << "Dynamic load balancing can only be activated for parallel "
            << "simulations." << nl
            << "Please run the simulation in parallel"
            << exit(FatalError);
       

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    // In the first iteration the cpuTimePerParticle_ is not yet set and
    // time statistics have to be gathered first
    if (firstTime_)
    {
        // First create the local cell list
        buildCellDataList(deltaT);

        solveCellList(cellDataList_);

        // Update the cell values
        forAll(cellDataList_,celli)
        {
            #ifdef FULLDEBUG
            if (cellDataList_[celli].cellID() != celli)
                FatalError << "cellID of cellDataList is "
                           << cellDataList_[celli].cellID() << " and does not "
                           << "match " << celli << exit(FatalError);
            #endif
            
            auto& cellData = cellDataList_[celli];

            this->deltaTChem_[celli] = cellData.deltaTChem();

            // Copy over the results
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<this->nSpecie_; i++)
            {
                this->RR_[i][celli] = cellData.RR()[i];
            }
        }
        
        updateTotalCpuTime(cellDataList_);

        firstTime_=false;

        return deltaTMin;
    }

    // If it is not the first time, the cellDataList has to be updated
    updateCellDataList(deltaT);

    // Get percentage of particles to send/receive from other processors
    auto sendAndReceiveData = getProcessorBalancing();

    const List<Tuple2<scalar,label>>& sendDataInfo = sendAndReceiveData.first();
    const List<label>& recvProc = sendAndReceiveData.second();

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    
    // indices of the reactCellList to create the sub lists to send 
    label start = 0;
    
    // Send all particles 
    for (auto& sendDataInfoI : sendDataInfo)
    {
        const scalar percToSend = sendDataInfoI.first();
        const label toProc = sendDataInfoI.second();
        
        label end = cellDataList_.size();

        cellsToSend
        (
            cellDataList_,
            totalCpuTime_*percToSend,
            start,
            end
        );
        
        // send particles 
        UOPstream toBuffer(toProc,pBufs);
        label dataSize = end - start;
        toBuffer << dataSize;
        for (label i=start; i < end; i++)
        {
            toBuffer << cellDataList_[i];
        }
        
        start = end;
    }
    
    // Set local to compute particle list
    SubList<baseDataContainer> localToComputeParticles
    (
        cellDataList_,
        cellDataList_.size()-start,
        start
    );

    pBufs.finishedSends();

    DynamicList<baseDataContainer> processorCells;
    
    List<label> receivedDataSizes(recvProc.size());
    
    // Read the received information
    forAll(recvProc,i)
    {
        label procI = recvProc[i];
        UIPstream fromBuffer(procI,pBufs);
        label dataSize;
        fromBuffer >> dataSize;
        receivedDataSizes[i] = dataSize;
        processorCells.reserve(processorCells.size()+dataSize);
        
        for (label k=0; k < dataSize; k++)
        {
            baseDataContainer p(fromBuffer);
            processorCells.append(std::move(p));
        }
    }    

    // Start solving local to compute particles
    solveCellList(localToComputeParticles);

    // Solve the chemistry on processor particles
    solveCellList(processorCells);

    // Send the information back 
    // Note: Now the processors to which we originally had send informations
    //       are the ones we receive from and vice versa 
    
    pBufs.clear();
    
    label pI = 0;
    
    forAll(recvProc,i)
    {
        label procI = recvProc[i];
        UOPstream toBuffer(procI,pBufs);
        
        toBuffer << receivedDataSizes[i];
        
        for (label k=0; k < receivedDataSizes[i]; k++)
            toBuffer << processorCells[pI++];
    }
    
    pBufs.finishedSends();
    
    start = 0;
    // Receive the particles --> now the sendDataInfo becomes the receive info
    for (auto& sendDataInfoI : sendDataInfo)
    {
        const label fromProc = sendDataInfoI.second();

        // send particles 
        UIPstream fromBuffer(fromProc,pBufs);
        label dataSize;
        fromBuffer >> dataSize;
        
        label end = start + dataSize;
        
        for (label i=start; i < end; i++)
            fromBuffer >> cellDataList_[i];
        
        start = end;
    }

    // Update the cell values
    forAll(cellDataList_,celli)
    {
        #ifdef FULLDEBUG
        if (cellDataList_[celli].cellID() != celli)
            FatalError << "cellID of cellDataList is "
                        << cellDataList_[celli].cellID() << " and does not "
                        << "match " << celli << exit(FatalError);
        #endif
            
        auto& cellData = cellDataList_[celli];

        this->deltaTChem_[celli] = cellData.deltaTChem();

        // Copy over the results
        deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

        this->deltaTChem_[celli] =
            min(this->deltaTChem_[celli], this->deltaTChemMax_);

        for (label i=0; i<this->nSpecie_; i++)
        {
            this->RR_[i][celli] = cellData.RR()[i];
        }
    }
    
    updateTotalCpuTime(cellDataList_); 

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}



