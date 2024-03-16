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

\*---------------------------------------------------------------------------*/

#include "TDACDataContainer.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * * 
Foam::TDACDataContainer::TDACDataContainer(Istream& is)
{
    read(is);
}

// * * * * * * * * * * * * * * * IO Functions * * * * * * * * * * * * * * * * *

Foam::Istream& Foam::TDACDataContainer::read
(
    Istream& is
)
{    
    phiq_.clear();
    c_.clear();
    c0_.clear();

    // Read scalar lists
    is >> phiq_;
    is >> c_;
    is >> c0_;

    // Read scalars
    is >> T_;
    is >> p_;
    is >> rho_;
    is >> deltaT_;
    is >> deltaTChem_;
    is >> cpuTime_;
    is >> addToTableCpuTime_;
    is >> procID_;
    is >> cellID_;

    is >> local_;

    return is;
}

Foam::Ostream& Foam::TDACDataContainer::write
(
    Ostream& os
) const
{
    // When the particle container is transferred or written it is no longer
    // local
    bool localParticle = false;

    os << phiq_ << token::SPACE
       << c_ << token::SPACE
       << c0_ << token::SPACE
       << T_ << token::SPACE
       << p_ << token::SPACE
       << rho_ << token::SPACE
       << deltaT_ << token::SPACE
       << deltaTChem_ << token::SPACE 
       << cpuTime_ << token::SPACE
       << addToTableCpuTime_ << token::SPACE
       << procID_ << token::SPACE
       << cellID_ << token::SPACE
       << localParticle << endl;
    return os;
}


Foam::Istream& Foam::operator >>(Istream& is, TDACDataContainer& eField)
{
    return eField.read(is);
}


Foam::Ostream& Foam::operator <<(Ostream& os, const TDACDataContainer& eField)
{
    return eField.write(os);
}
