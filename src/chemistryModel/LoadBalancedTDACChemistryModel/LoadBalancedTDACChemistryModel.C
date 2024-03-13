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

#include "LoadBalancedTDACChemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"



template<class ReactionThermo, class ThermoType>
Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::LoadBalancedTDACChemistryModel
(
    ReactionThermo& thermo
)
:
    TDACChemistryModel<ReactionThermo, ThermoType>(thermo)
{
    // Create the table for the remote cell computation
    tabulationRemote_ = chemistryTabulationMethod<ReactionThermo, ThermoType>::New
    (
        *this,
        *this
    );

    auto dict = this->subDictOrAdd("LoadBalancedTDACCoeffs");

    maxIterUpdate_ = dict.template getOrDefault<label>("updateIter",0);
    Info << "updateIter: "<<maxIterUpdate_<<endl;

    // Set iter to maxIterUpdate to force update in the first iteration
    iter_ = maxIterUpdate_;


    // Initialize cell list
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    cellDataField_.reserve(rho.size());

    // Number of additional properties stored in phiq.
    // Default is pressure and temperatur, if variableTimeStep is active
    // deltaT is added as well
    label nAdditions = 2;
    if (this->tabulation_->variableTimeStep())
        nAdditions = 3;

    forAll(rho,celli)
    {
        TDACDataContainer cData(this->nSpecie_,nAdditions);
        cellDataField_.append(cData);
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>::addCell
(
    DynamicList<TDACDataContainer*>& cellList,
    const scalarField& phiq,
    const scalar& T,
    const scalar& p,
    const scalar& rho,
    const scalar& deltaT,
    const label& celli
)
{
    label MyProcNo = Pstream::myProcNo();

    // Checkout the cData container from the field
    TDACDataContainer& cData = cellDataField_[celli];

    cData.phiq() = phiq;

    cData.T() = T;

    cData.p() = p;

    cData.rho() = rho;

    cData.deltaT() = deltaT;

    cData.deltaTChem() = this->deltaTChem_[celli];

    cData.proc() = MyProcNo;

    cData.cellID() = celli;

    cellList.append(&cData);

    cellsToSolve_++;
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::cellsToSend
(
    const DynamicList<TDACDataContainer*>& cellList,
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

        cpuTime += cellList[i]->cpuTime()+cellList[i]->addToTableCpuTime();
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::updateProcessorBalancing()
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

    // list of the distributed load for all processors
    List<DynamicList<sendDataStruct>> distributedLoadAllProcs(numProcs);
    
    // list of processors to receive data from
    List<DynamicList<label>> receiveDataFromProc(numProcs);

    // balance the load by calculating the percentages to be send 
    forAll(sortedCpuTimeOnProcessors,i)
    {
        const label procI = sortedCpuTimeOnProcessors[i].second[0];

        // List of processors to send information to
        DynamicList<sendDataStruct>& sendLoadList = 
            distributedLoadAllProcs[procI];
        
        // Reserve space
        sendLoadList.reserve
        (
            std::floor(0.5*(numProcs-i))
        );


        // Total time cpu time on procI -- including time to solve ODE and 
        // adding to the table
        const scalar cpuTimeProcI = sortedCpuTimeOnProcessors[i].first;
        
        scalar cpuTimeOverhead = cpuTimeProcI - averageCpuTime;

        if (cpuTimeOverhead < 0)
            continue;

        // Loop over the other processors and distribute the load
        // in reverse order
        for (label k=numProcs-1; k > 0; k--)
        {
            const label procK = sortedCpuTimeOnProcessors[k].second[0];
            
            const scalar cpuTimeProcK = sortedCpuTimeOnProcessors[k].first;
            
            const scalar capacityOfProcK = averageCpuTime - cpuTimeProcK;
            
            if (capacityOfProcK <= 0)
                continue;
            
            const scalar newCapacity = capacityOfProcK - cpuTimeOverhead;
            
            // Only send information to procK if it is larger than 2% of the
            // total average cell time
            if ((std::min(capacityOfProcK,cpuTimeOverhead)/cpuTimeProcI) < 0.02)
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
                    sendDataStruct
                    (
                        cpuTimeOverhead,
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
                    sendDataStruct
                    (
                        capacityOfProcK,
                        sortedCpuTimeOnProcessors[k].second[0]
                    )
                );
            }
        }
    }

    sendAndReceiveData_.first() = distributedLoadAllProcs[Pstream::myProcNo()];
    sendAndReceiveData_.second() = receiveDataFromProc[Pstream::myProcNo()];

    // Create send and receive lists
    sendToProcessor_.resize(Pstream::nProcs());
    receiveFromProcessor_.resize(Pstream::nProcs());

    // Set all to false
    forAll(sendToProcessor_,i)
    {
        sendToProcessor_[i] = false;
        receiveFromProcessor_[i] = false;
    }

    for (auto& sendInfoData : sendAndReceiveData_.first())
    {
        sendToProcessor_[sendInfoData.toProc] = true;
    }

    for (label procID : sendAndReceiveData_.second())
    {
        receiveFromProcessor_[procID] = true;
    }
}


template<class ReactionThermo, class ThermoType>
Foam::List<std::pair<scalar,Foam::List<scalar>>> 
Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::getSortedCPUTimesOnProcessor() const
{
    // Number of processors 
    const scalar numProcs = Pstream::nProcs();
        
    // Gather the data from all processors
    List<List<scalar>> cpuTimeOnProcessors(numProcs);
    cpuTimeOnProcessors[Pstream::myProcNo()].resize(4);
    cpuTimeOnProcessors[Pstream::myProcNo()][0] = totalCpuTime_;
    cpuTimeOnProcessors[Pstream::myProcNo()][1] = cellsToSolve_;
    cpuTimeOnProcessors[Pstream::myProcNo()][2] = addToTableCpuTime_;
    cpuTimeOnProcessors[Pstream::myProcNo()][3] = searchISATCpuTime_;


    Pstream::gatherList(cpuTimeOnProcessors);
    Pstream::scatterList(cpuTimeOnProcessors);

    // use std::pair for std::sort algorithm
    // Data structure is:
    // List[index].first:  CPU time required for last time step
    // List[index].second: Returns a List<scalar> with:
    //                       List[0]: processor ID
    //                       List[1]: Number of cells to compute 
    //                       List[2]: Time to add cells to table
    List<std::pair<scalar,List<scalar>>> sortedCpuTimeOnProcessors(numProcs);
    
    forAll(sortedCpuTimeOnProcessors,i)
    {
        sortedCpuTimeOnProcessors[i].first = cpuTimeOnProcessors[i][0];
        sortedCpuTimeOnProcessors[i].second.resize(4);
        sortedCpuTimeOnProcessors[i].second[0] = i;
        sortedCpuTimeOnProcessors[i].second[1] = cpuTimeOnProcessors[i][1];
        sortedCpuTimeOnProcessors[i].second[2] = cpuTimeOnProcessors[i][2];
        sortedCpuTimeOnProcessors[i].second[3] = cpuTimeOnProcessors[i][3];
    }
    
    // sort using std::sort()
    // first entry has largest cpu time --> descending order
    std::stable_sort
    (
        sortedCpuTimeOnProcessors.begin(), 
        sortedCpuTimeOnProcessors.end(),
        std::greater<std::pair<scalar,List<scalar>>>()
    );

    return sortedCpuTimeOnProcessors;
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::updateTotalCpuTime
(
    const DynamicList<TDACDataContainer*>& cellList
)
{
    // Calculate the total time spent solving the particles on this processor
    totalCpuTime_  = 0;
    addToTableCpuTime_ = 0;
    for (const auto& cDataPtr : cellList)
    {
        addToTableCpuTime_ += cDataPtr->addToTableCpuTime();
        totalCpuTime_ += (cDataPtr->cpuTime()+cDataPtr->addToTableCpuTime());
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>::solveCell
(
    TDACDataContainer& cData
)
{
    // Store total time waiting to attribute to add or grow
    scalar timeTmp = clockTime_.timeIncrement();

    // Note: first nSpecie entries are the Yi values in phiq
    Field<scalar>& c = cData.c();
    Field<scalar>& c0 = cData.c0();
    for (label i=0; i<this->nSpecie_; i++)
    {
        c[i] = cData.rho()*cData.phiq()[i]/this->specieThermo_[i].W();
        c0[i] = c[i];
    }

    const bool reduced = this->mechRed()->active();

    scalar timeLeft = cData.deltaT();

    scalar reduceMechCpuTime_ = 0;


    if (reduced)
    {
        // Reduce mechanism change the number of species (only active)
        this->mechRed()->reduceMechanism(c, cData.T(), cData.p());
        scalar timeIncr = clockTime_.timeIncrement();
        reduceMechCpuTime_ += timeIncr;
        timeTmp += timeIncr;
        nSpecieReduced_ = this->nSpecie_;
    }

    // Calculate the chemical source terms
    while (timeLeft > SMALL)
    {
        scalar dt = timeLeft;
        if (reduced)
        {
            // completeC_ used in the overridden ODE methods
            // to update only the active species
            this->completeC_ = c;

            // Solve the reduced set of ODE
            this->solve
            (
                this->simplifiedC_, cData.T(), cData.T(), dt, cData.deltaTChem()
            );

            for (label i=0; i<this->NsDAC_; ++i)
            {
                c[this->simplifiedToCompleteIndex_[i]] = this->simplifiedC_[i];
            }
        }
        else
        {
            this->solve(c, cData.T(), cData.p(), dt, cData.deltaTChem());
        }
        timeLeft -= dt;
    }

    {
        scalar timeIncr = clockTime_.timeIncrement();
        solveChemistryCpuTime_ += timeIncr;
    }

    // When operations are done and if mechanism reduction is active,
    // the number of species (which also affects nEqns) is set back
    // to the total number of species (stored in the this->mechRed object)
    if (reduced)
    {
        this->nSpecie_ = this->mechRed()->nSpecie();
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::solveCellList
(
    UList<TDACDataContainer*>& cellList,
    const bool isLocal
)
{
    for (TDACDataContainer* cDataPtr : cellList)
    {
        // Check if it now can be found in table
        // this is possible if a previous cell computed this result
        if (lookUpCellInTable(*cDataPtr,isLocal))
            continue;

        // We cannot use here cpuTimeIncrement() of OpenFOAM as this 
        // returns only measurements in 100Hz or 1000Hz intervals depending
        // on the installed kernel 
        auto start = std::chrono::high_resolution_clock::now();

        solveCell
        (
            *cDataPtr
        );

        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end-start));

        // Add the time required to solve this particle to the list 
        // as seconds
        cDataPtr->cpuTime() = duration.count()*1.0E-6;

        // Add to table
        // Does not recompute the reduced reaction mechanism as it was just 
        // computed for this cell in solveCell()
        addCellToTable(*cDataPtr,isLocal,false);
    }
}


template<class ReactionThermo, class ThermoType>
bool Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::lookUpCellInTable
(
    TDACDataContainer& cData,
    const bool isLocal
)
{
    chemistryTabulationMethod<ReactionThermo, ThermoType>* tabPtr;
    if (isLocal)
        tabPtr = this->tabulation_.get();
    else
        tabPtr = this->tabulationRemote_.get();

    if 
    (
        tabPtr->active() 
     && tabPtr->retrieve(cData.phiq(), Rphiq_)
    )
    {
        // Note: first nSpecie entries are the Yi values in phiq
        Field<scalar>& c = cData.c();
        Field<scalar>& c0 = cData.c0();
        for (label i=0; i<this->nSpecie_; i++)
        {
            c0[i] = cData.rho()*cData.phiq()[i]/this->specieThermo_[i].W();
        }

        // Retrieved solution stored in Rphiq_
        for (label i=0; i<this->nSpecie(); ++i)
        {
            c[i] = cData.rho()*Rphiq_[i]/this->specieThermo_[i].W();
        }
        return true;
    }

    return false;
}


template<class ReactionThermo, class ThermoType>
void Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::addCellToTable
(
    const TDACDataContainer& cData,
    const bool isLocal,
    const bool requiresRecomputeReducedMech
)
{
    chemistryTabulationMethod<ReactionThermo, ThermoType>* tabPtr;

    if (isLocal)
        tabPtr = this->tabulation_.get();
    else
        tabPtr = this->tabulationRemote_.get();

    // We cannot use here cpuTimeIncrement() of OpenFOAM as this 
    // returns only measurements in 100Hz or 1000Hz intervals depending
    // on the installed kernel 
    auto start = std::chrono::high_resolution_clock::now();

    const auto& c = cData.c();

    // Not sure if this is necessary
    Rphiq_ = Zero;

    forAll(c, i)
    {
        Rphiq_[i] = c[i]/cData.rho()*this->specieThermo_[i].W();
    }
    if (tabPtr->variableTimeStep())
    {
        Rphiq_[Rphiq_.size()-3] = cData.T();
        Rphiq_[Rphiq_.size()-2] = cData.p();
        Rphiq_[Rphiq_.size()-1] = cData.deltaT();
    }
    else
    {
        Rphiq_[Rphiq_.size()-2] = cData.T();
        Rphiq_[Rphiq_.size()-1] = cData.p();
    }

    clockTime_.timeIncrement();

    // If tabulation is used, we add the information computed here to
    // the stored points (either expand or add)
    if 
    (
        tabPtr->active() 
        && !tabPtr->retrieve(cData.phiq(), Rphiq_)
    )
    {
        if (this->mechRed()->active())
        {
            if (requiresRecomputeReducedMech)
                this->mechRed_->reduceMechanism(c, cData.T(), cData.p());
            else
                this->setNSpecie(nSpecieReduced_);
        }

        label growOrAdd =
            tabPtr->add
            (
                cData.phiq(), Rphiq_, cData.rho(), cData.deltaT()
            );

        // Only collect information for local cells
        if (isLocal)
        {
            const label celli = cData.cellID();

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end-start));
            // Add the time to the cell data container
            cData.setAddToTableCpuTime(duration.count()*1.0E-6);

            if (growOrAdd)
            {
                this->setTabulationResultsAdd(celli);
                addNewLeafCpuTime_ += clockTime_.timeIncrement();
            }
            else
            {
                this->setTabulationResultsGrow(celli);
                growCpuTime_ += clockTime_.timeIncrement();
            }
        }
    
        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the this->mechRed object)
        if (this->mechRed()->active())
        {
            this->nSpecie_ = this->mechRed()->nSpecie();
        }
    }
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::solve
(
    const DeltaTType& deltaT
)
{
    // List of cells that need to be solved
    DynamicList<TDACDataContainer*> cellList(cellDataField_.size());

    cellsToSolve_ = 0;

    // Increment counter of time-step
    this->timeSteps_++;

    const bool reduced = this->mechRed()->active();

    label nAdditionalEqn = (this->tabulation_->variableTimeStep() ? 1 : 0);

    basicSpecieMixture& composition = this->thermo().composition();

    // CPU time analysis
    clockTime_= clockTime();
    clockTime_.timeIncrement();
    reduceMechCpuTime_ = 0;
    addNewLeafCpuTime_ = 0;
    growCpuTime_ = 0;
    solveChemistryCpuTime_ = 0;
    searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    Rphiq_.resize(this->nEqns() + nAdditionalEqn);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<this->nSpecie_; i++)
        {
            const volScalarField& Yi = this->Y_[i];
            c[i] = rhoi*Yi[celli]/this->specieThermo_[i].W();
            c0[i] = c[i];
            phiq[i] = Yi[celli];
        }
        phiq[this->nSpecie()] = Ti;
        phiq[this->nSpecie() + 1] = pi;
        if (this->tabulation_->variableTimeStep())
        {
            phiq[this->nSpecie() + 2] = deltaT[celli];
        }

        // Not sure if this is necessary
        Rphiq_ = Zero;

        clockTime_.timeIncrement();

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method
        if (this->tabulation_->active() && this->tabulation_->retrieve(phiq, Rphiq_))
        {
            // Retrieved solution stored in Rphiq_
            for (label i=0; i<this->nSpecie(); ++i)
            {
                c[i] = rhoi*Rphiq_[i]/this->specieThermo_[i].W();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        
            // Set the RR vector (used in the solver)
            for (label i=0; i<this->nSpecie_; ++i)
            {
                this->RR_[i][celli] =
                    (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
            }
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information for later compuation
        // with the load balanced method
        else
        {
            addCell(cellList,phiq,Ti,pi,rhoi,deltaT[celli],celli);
        }
    }

    // Update total cpu time based on the current cell list
    updateTotalCpuTime(cellList); 

    // ========================================================================
    // Solve for cells that were not found in the table
    // ========================================================================

    // If it is solved the first time, computational statistics have to be 
    // gathered first
    if (firstTime_)
    {
        solveCellList(cellList,true);
        firstTime_ = false;
    }
    else
    {
        if (iter_++ >= maxIterUpdate_)
        {
            updateProcessorBalancing();
            iter_ = 0;
        }

        const List<sendDataStruct>& sendDataInfo = sendAndReceiveData_.first();
        const List<label>& recvProc = sendAndReceiveData_.second();
        
        // indices of the reactCellList to create the sub lists to send 
        label start = 0;
        
        // Send all particles 
        for (auto& sendDataInfoI : sendDataInfo)
        {
            const scalar cpuTimeToSend = sendDataInfoI.cpuTimeToSend;
            const label toProc = sendDataInfoI.toProc;
            
            label end = cellList.size();

            cellsToSend
            (
                cellList,
                cpuTimeToSend,
                start,
                end
            );
            
            // send particles 
            UOPstream toBuffer
            (
                pBufs_.commsType(),
                toProc,
                pBufs_.sendBuffer(toProc),
                pBufs_.tag(),
                pBufs_.comm(),
                false
            );

            label dataSize = end - start;
            toBuffer << dataSize;
            for (label i=start; i < end; i++)
            {
                toBuffer << *(cellList[i]);
            }
            
            start = end;
        }
        
        // Set local to compute particle list
        SubList<TDACDataContainer*> localToComputeCells
        (
            cellList,
            cellList.size()-start,
            start
        );

        pBufs_.finishedSends(sendToProcessor_,receiveFromProcessor_);

        DynamicList<TDACDataContainer> processorCells;
        
        List<label> receivedDataSizes(recvProc.size());
        
        // Read the received information
        forAll(recvProc,i)
        {
            label receiveBufferPosition=0;
            label procI = recvProc[i];
            UIPstream fromBuffer
            (
                pBufs_.commsType(),
                procI,
                pBufs_.recvBuffer(procI),
                receiveBufferPosition,
                pBufs_.tag(),
                pBufs_.comm(),
                false
            );

            label dataSize;
            fromBuffer >> dataSize;
            receivedDataSizes[i] = dataSize;
            processorCells.reserve(processorCells.size()+dataSize);
            
            for (label k=0; k < dataSize; k++)
            {
                TDACDataContainer p(fromBuffer);
                processorCells.append(std::move(p));
            }
        }    

        // We need to create a pointer list for the solveCellList function
        List<TDACDataContainer*> processorCellsPtr(processorCells.size(), nullptr);

        // We cannot create this pointer list in the for loop before, as each
        // resize of the dynamic list will invalidate the pointers. 
        forAll(processorCells,i)
        {
            processorCellsPtr[i] = &processorCells[i];
        }

        // Solve the local cells first
        solveCellList(localToComputeCells,true);

        // Solve the chemistry on processor particles
        solveCellList(processorCellsPtr,false);

        // Send the information back 
        // Note: Now the processors to which we originally had send informations
        //       are the ones we receive from and vice versa 
        
        pBufs_.switchSendRecv();
        
        label pI = 0;
        
        forAll(recvProc,i)
        {
            label procI = recvProc[i];
            UOPstream toBuffer
            (
                pBufs_.commsType(),
                procI,
                pBufs_.sendBuffer(procI),
                pBufs_.tag(),
                pBufs_.comm(),
                false
            );
            
            toBuffer << receivedDataSizes[i];
            
            for (label k=0; k < receivedDataSizes[i]; k++)
                toBuffer << processorCells[pI++];
        }
        
        pBufs_.finishedSends();
        
        start = 0;
        label end = 0;
        // Receive the particles --> now the sendDataInfo becomes the receive info
        for (auto& sendDataInfoI : sendDataInfo)
        {
            const label fromProc = sendDataInfoI.toProc;

            // send particles 
            label receiveBufferPosition=0;
            label dataSize;
            UIPstream fromBuffer
            (
                pBufs_.commsType(),
                fromProc,
                pBufs_.recvBuffer(fromProc),
                receiveBufferPosition,
                pBufs_.tag(),
                pBufs_.comm(),
                false
            );

            fromBuffer >> dataSize;
            
            end = start + dataSize;
            
            for (label i=start; i < end; i++)
                fromBuffer >> *(cellList[i]);
            
            start = end;
        }

        // Switch sendAndRecv back
        pBufs_.switchSendRecv();

        // =====================================================================
        //                      Update Table for Remote cells
        // =====================================================================

        // Remote cells are all from start to end 
        SubList<TDACDataContainer*> remoteComputedCells
        (
            cellList,
            0,
            end
        );

        for (const auto& cDataPtr : remoteComputedCells)
        {
            // dereference pointer
            const auto& cData = *cDataPtr;

            // Add cell to ISAT table and log CPU time
            addCellToTable(cData,true);
        }
    }

    // ========================================================================
    //                      Update Reaction Rate
    // ========================================================================

    for (const auto& cDataPtr : cellList)
    {
        // dereference pointer
        const auto& cData = *cDataPtr;

        const label celli = cData.cellID();

        const auto& c = cData.c();
        const auto& c0 = cData.c0();

        deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

        this->deltaTChem_[celli] =
            min(this->deltaTChem_[celli], this->deltaTChemMax_);

        // Set the RR vector (used in the solver)
        for (label i=0; i<this->nSpecie_; ++i)
        {
            this->RR_[i][celli] =
                (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
        }
    }

    if (this->mechRed_->log() || this->tabulation_->log())
    {
        this->cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (this->mechRed_->log())
    {
        this->cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (this->tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        this->tabulation_->update();
        this->tabulationRemote_->update();

        // Write the performance of the tabulation
        this->tabulation_->writePerformance();

        if (this->tabulation_->log())
        {
            this->cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            this->cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            this->cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && this->mechRed_->log())
    {
        // Write average number of species
        this->nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies/nAvg << endl;
    }

    if (reduced && Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineReduce(active, orEqOp<bool>());

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    forAll(this->Y(), i)
    {
        if (composition.active(i))
        {
            this->Y()[i].writeOpt(IOobject::AUTO_WRITE);
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::solve
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
Foam::scalar Foam::LoadBalancedTDACChemistryModel<ReactionThermo, ThermoType>
::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}

// ************************************************************************* //
