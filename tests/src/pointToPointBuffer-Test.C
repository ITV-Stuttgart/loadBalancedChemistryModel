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
    Test the pointToPointBuffer 

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
#include "pointToPointBuffer.H"


TEST_CASE("pointToPointBuffer-Test","[Pstream]")
{
    // =========================================================================
    //                      Prepare Case
    // =========================================================================
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    #include "createMesh.H"

    Info << "Create buffer"<<endl;
    // Create a pointToPointBuffer
    pointToPointBuffer pBuf;

    // Test matrix requires execution with 4 cores!
    // Processor 0 sends a list to processor 3 with 10 label entries
    // And Processor 3 sends processor 0 a list with 5 label entries


    labelList myData;
    if (Pstream::myProcNo() == 0)
    {
        myData.resize(10);
        forAll(myData,i)
        {
            myData[i] = 10;
        }
    }
    else if  (Pstream::myProcNo() == 3)
    {
        myData.resize(5);
        forAll(myData,i)
        {
            myData[i] = 5;
        }
    }


    // Create the UOPstream object
    if (Pstream::myProcNo() == 0)
    {
        label toProc = 3;
        UOPstream toBuffer(pBuf.commsType(),toProc,pBuf.sendBuffer(toProc),pBuf.tag(),pBuf.comm(),false);

        toBuffer << myData;
    }
    else if (Pstream::myProcNo() == 3)
    {
        label toProc = 0;
        UOPstream toBuffer(pBuf.commsType(),toProc,pBuf.sendBuffer(toProc),pBuf.tag(),pBuf.comm(),false);

        toBuffer << myData;
    }

    // Update the send and receive buffers
    pBuf.update();

    pBuf.finishedSends();

    // Now read data 
    // Create the UOPstream object
    labelList recvData;
    if (Pstream::myProcNo() == 0)
    {
        label fromProc = 3;
        label receiveBufferPosition=0;
        UIPstream fromBuffer(pBuf.commsType(),fromProc,pBuf.recvBuffer(fromProc),receiveBufferPosition,pBuf.tag(),pBuf.comm(),false);

        fromBuffer >> recvData;
        REQUIRE(recvData.size() == 5);
        forAll(recvData,i)
        {
            REQUIRE(recvData[i] == 5);
        }
    }
    else if (Pstream::myProcNo() == 3)
    {
        label fromProc = 0;
        label receiveBufferPosition=0;
        UIPstream fromBuffer(pBuf.commsType(),fromProc,pBuf.recvBuffer(fromProc),receiveBufferPosition,pBuf.tag(),pBuf.comm(),false);

        fromBuffer >> recvData;
        REQUIRE(recvData.size() == 10);
        forAll(recvData,i)
        {
            REQUIRE(recvData[i] == 10);
        }
    }

    // =========================================================================
    // Now send again from processor 0 to 3 and 3 to 0
    // But this time use double the number of entries
    Info << "*** Test with exchange sizes only ***"<<endl;


    if (Pstream::myProcNo() == 0)
    {
        myData.resize(15);
        forAll(myData,i)
        {
            myData[i] = 15;
        }
    }
    else if  (Pstream::myProcNo() == 3)
    {
        myData.resize(8);
        forAll(myData,i)
        {
            myData[i] = 8;
        }
    }

    // Use send/recv Lists
    List<bool> sendToProcessor(Pstream::nProcs(),false);
    List<bool> receiveFromProcessor(Pstream::nProcs(),false);

    // Create the UOPstream object
    if (Pstream::myProcNo() == 0)
    {
        label toProc = 3;
        UOPstream toBuffer(pBuf.commsType(),toProc,pBuf.sendBuffer(toProc),pBuf.tag(),pBuf.comm(),false);

        toBuffer << myData;
        sendToProcessor[toProc] = true;
        // Here we receive from the processor we send to
        receiveFromProcessor[toProc] = true; 
    }
    else if (Pstream::myProcNo() == 3)
    {
        label toProc = 0;
        UOPstream toBuffer(pBuf.commsType(),toProc,pBuf.sendBuffer(toProc),pBuf.tag(),pBuf.comm(),false);

        toBuffer << myData;
        sendToProcessor[toProc] = true;
        // Here we receive from the processor we send to
        receiveFromProcessor[toProc] = true; 
    }

    pBuf.finishedSends(sendToProcessor,receiveFromProcessor);

    // Now read data 
    // Create the UOPstream object
    recvData.clear();
    if (Pstream::myProcNo() == 0)
    {
        label fromProc = 3;
        label receiveBufferPosition=0;
        UIPstream fromBuffer(pBuf.commsType(),fromProc,pBuf.recvBuffer(fromProc),receiveBufferPosition,pBuf.tag(),pBuf.comm(),false);

        fromBuffer >> recvData;
        REQUIRE(recvData.size() == 8);
        forAll(recvData,i)
        {
            REQUIRE(recvData[i] == 8);
        }
    }
    else if (Pstream::myProcNo() == 3)
    {
        label fromProc = 0;
        label receiveBufferPosition=0;
        UIPstream fromBuffer(pBuf.commsType(),fromProc,pBuf.recvBuffer(fromProc),receiveBufferPosition,pBuf.tag(),pBuf.comm(),false);

        fromBuffer >> recvData;
        REQUIRE(recvData.size() == 15);
        forAll(recvData,i)
        {
            REQUIRE(recvData[i] == 15);
        }
    }
}

