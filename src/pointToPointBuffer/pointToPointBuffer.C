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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "pointToPointBuffer.H"


void Foam::pointToPointBuffer::update(const List<label>& sendSizes)
{
    sendBufferSize_ = sendSizes;

    if (sendBufferSize_.size() != UPstream::nProcs(comm_))
    {
        FatalErrorInFunction
            << "Size of container " << sendBufferSize_.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm_)
            << Foam::abort(FatalError);
    }
 

    recvBufferSize_.resize_nocopy(sendBufferSize_.size());

    // Use UPstream::allToAll as it is used in exchangeSizes()
    UPstream::allToAll(sendBufferSize_, recvBufferSize_, comm_);
    
    // Update the receive buffer sizes
    forAll(recvBufferList_,procI)
    {
        recvBufferList_[procI].resize(recvBufferSize_[procI]);
    }
}


void Foam::pointToPointBuffer::update()
{
    forAll(sendBufferList_,procI)
    {
        sendBufferSize_[procI] = sendBufferList_[procI].byteSize();
    }

    recvBufferSize_.resize_nocopy(sendBufferSize_.size());

    // Use UPstream::allToAll as it is used in exchangeSizes()
    UPstream::allToAll(sendBufferSize_, recvBufferSize_, comm_);
    
    // Update the receive buffer sizes
    forAll(recvBufferList_,procI)
    {
        recvBufferList_[procI].resize(recvBufferSize_[procI]);
    }
}


void Foam::pointToPointBuffer::switchSendRecv()
{
    labelList recvBufferSizeCopy = recvBufferSize_;
    forAll(sendBufferSize_,procI)
    {
        recvBufferSize_[procI] = sendBufferSize_[procI];
        sendBufferSize_[procI] = recvBufferSizeCopy[procI];
    }
    // Update the receive buffer sizes
    forAll(recvBufferList_,procI)
    {
        recvBufferList_[procI].resize(recvBufferSize_[procI]);
    }
}


void Foam::pointToPointBuffer::finishedSends()
{
    const label startOfRequests = UPstream::nRequests();

    checkBufferSize();

    forAll(sendBufferSize_,procI)
    {
        if (recvBufferSize_[procI] > 0 && procI != Pstream::myProcNo())
        {
            IPstream::read
            (
                commsType_,
                procI,
                recvBufferList_[procI].data(),
                recvBufferSize_[procI],
                UPstream::msgType(),
                comm_
            );  
        }
    }
    
    forAll(sendBufferSize_,procI)
    {
        if (sendBufferSize_[procI] > 0 && procI != Pstream::myProcNo())
        {
            OPstream::write
            (
                commsType_,
                procI,
                sendBufferList_[procI].cdata(),
                sendBufferList_[procI].byteSize(),
                UPstream::msgType(),
                comm_
            );
        }
    }

    UPstream::waitRequests(startOfRequests);

    // Clear the send buffer
    forAll(sendBufferList_,procI)
    {
        sendBufferList_[procI].clear();
    }
}


void Foam::pointToPointBuffer::checkBufferSize()
{
    forAll(sendBufferSize_,procI)
    {
        if (sendBufferList_[procI].byteSize() != sendBufferSize_[procI])
            FatalError 
                << "sendBufferSize "<<sendBufferList_[procI].byteSize()
                << " for processor "<<procI 
                << " does not match " << sendBufferSize_[procI] << ". Consider "
                << "calling update() to set the new buffer sizes" 
                << exit(FatalError);

        if (recvBufferList_[procI].byteSize() != recvBufferSize_[procI])
            FatalError 
                << "recvBufferSize "<<recvBufferList_[procI].byteSize()
                << " for processor "<<procI 
                << " does not match " << recvBufferSize_[procI] << ". Consider "
                << "calling update() to set the new buffer sizes" 
                << exit(FatalError);
    }
}

