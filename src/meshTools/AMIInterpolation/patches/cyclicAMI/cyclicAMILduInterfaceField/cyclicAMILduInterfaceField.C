/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "cyclicAMILduInterfaceField.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cyclicAMILduInterfaceField, 0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMILduInterfaceField::~cyclicAMILduInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicAMILduInterfaceField::transformCoupleField
(
    solveScalarField& f,
    const direction cmpt
) const
{
    if (doTransform())
    {
        foamExecutor exec;
        auto f_p = f.begin();
        const auto forwardT_p = forwardT().cbegin();
        const auto localRank = rank();
        if (forwardT().size() == 1)
        {
            auto Lambda = [=](label i){
                f_p[i] *=  pow(diag(forwardT_p[0]).component(cmpt), localRank);
            };
            exec.parallelFor(Lambda, f.size());
        }
        else
        {
            auto Lambda = [=](label i){
                f_p[i] *= pow(diag(forwardT_p[i]).component(cmpt), localRank);
            };
            exec.parallelFor(Lambda,f.size());
        }
    }
}


// ************************************************************************* //
