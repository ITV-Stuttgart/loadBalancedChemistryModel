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
    Foam::pointToPointBuffer
 
Description
    Provides a buffer for inter-processor communication which only exchanges
    sizes if the update() function is called. Otherwise the buffer sizes are 
    fixed, and no MPI_allToAll is required. In the OpenFOAM PstreamBuffers the 
    sizes are exchanged for each finishedSends which requires an allToAll 
    communication. Especially for large scale operations this is expansive.

    This is closer to the processorLduInterface class than PstreamBuffers.

SourceFiles
    pointToPointBuffer.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024

\*---------------------------------------------------------------------------*/

#ifndef pointToPointBuffer_H
#define pointToPointBuffer_H

#include "fvCFD.H"

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class pointToPointBuffer
\*---------------------------------------------------------------------------*/

class pointToPointBuffer
{
    protected:

    // Protected Member Variables

        //- Send buffer list
        //  Size of list is Pstream::nProcs()
        List<DynamicList<char>> sendBufferList_;

        //- Receive buffer list
        List<DynamicList<char>> recvBufferList_;
    
        //- Store the size of the buffer
        List<std::streamsize> sendBufferSize_;

        List<std::streamsize> recvBufferSize_;

        //- Always non-blocking communication
        const UPstream::commsTypes commsType_;

        //- World communicator
        const label comm_;


    // Protected Member Functions

        //- Check that the buffer is the expected size
        void checkBufferSize();

        //- Echange the send/recv buffer sizes based on the given send 
        //  and receive lists
        void exchangeBufferSizes
        (
            const List<bool>& sendToProcessor,
            const List<bool>& receiveFromProcessor
        );

    public:
    
        pointToPointBuffer()
        : 
            sendBufferList_(Pstream::nProcs()),
            recvBufferList_(Pstream::nProcs()),
            sendBufferSize_(Pstream::nProcs()),
            recvBufferSize_(Pstream::nProcs()),
            commsType_(Pstream::commsTypes::nonBlocking),
            comm_(UPstream::worldComm)
        {};   

    // Modify

        //- Update the send and receive buffers based on the current 
        //  sendBufferList
        void update();

        //- Exchange all information
        //  Requires update() to be called if the buffer changed
        void finishedSends();

        //- Send the buffers to the processors defined by sendToProcessor
        //  and receiveFromProcessor
        //  Uses exchangeBufferSizes
        void finishedSends
        (
            const List<bool>& sendToProcessor,
            const List<bool>& receiveFromProcessor
        );

        //- Switch send and receive buffer sizes
        //  Helps for cases when data is send for computation and then received
        //  again
        void switchSendRecv();


    // Access

        //- Return the recv buffer for procI
        DynamicList<char>& recvBuffer(const label procI)
        {return recvBufferList_[procI];}

        // Return the send buffer for procI
        DynamicList<char>& sendBuffer(const label procI)
        {return sendBufferList_[procI];}

        //- Return the communication type
        const UPstream::commsTypes& commsType() {return commsType_;}

        //- Return the tag
        const int& tag() {return UPstream::msgType();}

        //- Return communicator world
        const label& comm() {return comm_;}


};

}   // End of namespace Foam
#endif
