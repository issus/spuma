/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::flipOp

Description
    Functor to negate primitives. Dummy for most other types.

    Used in mesh transformations where face can flip.

SourceFiles
    flipOp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flipOp_H
#define Foam_flipOp_H

#include "fieldTypes.H"
#include "deviceM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noOp Declaration
\*---------------------------------------------------------------------------*/

//- Same as Foam::identityOp. Should \em never be specialized.
struct noOp
{
    template<class T>
    FOAM_DEVICE const T& operator()(const T& val) const noexcept
    {
        return val;
    }
};


/*---------------------------------------------------------------------------*\
                           Class flipOp Declaration
\*---------------------------------------------------------------------------*/

struct flipOp
{
    template<typename T>
    FOAM_DEVICE T operator()(const T& val) const
    {
        return val;
    }
};


//- Invert boolean value
struct flipBoolOp
{
    FOAM_DEVICE bool operator()(const bool val) const noexcept
    {
        return !val;
    }
};


//- Negate integer values
struct flipLabelOp
{
    FOAM_DEVICE label operator()(const label val) const
    {
        return -val;
    }
};


// Template specialisation for primitives that support negation
// template<> FOAM_DEVICE scalar flipOp::operator()(const scalar&) const;
// template<> FOAM_DEVICE vector flipOp::operator()(const vector&) const;
// template<> FOAM_DEVICE sphericalTensor flipOp::operator()(const sphericalTensor&) const;
// template<> FOAM_DEVICE symmTensor flipOp::operator()(const symmTensor&) const;
// template<> FOAM_DEVICE tensor flipOp::operator()(const tensor&) const;
// template<> FOAM_DEVICE triad flipOp::operator()(const triad&) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "flipOp.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
