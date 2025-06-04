/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "profiling.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolation::weightedSum
(
    const scalar lowWeightCorrection,
    const labelListList& allSlots,
    const scalarListList& allWeights,
    const scalarField& weightsSum,
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
)
{
    if (lowWeightCorrection > 0)
    {
        forAll(result, facei)
        {
            if (weightsSum[facei] < lowWeightCorrection)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& slots = allSlots[facei];
                const scalarList& weights = allWeights[facei];

                forAll(slots, i)
                {
                    cop(result[facei], facei, fld[slots[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& slots = allSlots[facei];
            const scalarList& weights = allWeights[facei];

            forAll(slots, i)
            {
                cop(result[facei], facei, fld[slots[i]], weights[i]);
            }
        }
    }
}

template<class Type, class CombineOp>
void Foam::AMIInterpolation::weightedSum
(
    const scalar lowWeightCorrection,
    const labelListList& allSlots,
    const scalarListList& allWeights,
    const labelList& flatSlots,
    const scalarList& flatWeights,
    const labelList& flatIdx,
    const scalarField& weightsSum,
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
)
{
//    List<Type> gpu_result(result);

/*
    //check whether flat and listlist are equal
    bool equal = true;
    forAll(allSlots,i)
    {
        const labelList& slots = allSlots[i];
        const scalarList& weights = allWeights[i];
        forAll(slots, j)
        {
            if
            (
                (slots[j] != flatSlots[flatIdx[i]+j]) ||
                (abs(weights[j]-flatWeights[flatIdx[i]+j]) > 1e-6)
            )
            {
                equal = false;
                break;
            }
        }

        if (!equal)
            break;

    }
    Info<<"flatIdx and ListList equal:"<< equal <<endl;

    if (lowWeightCorrection > 0)
    {
        forAll(result, facei)
        {
            if (weightsSum[facei] < lowWeightCorrection)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& slots = allSlots[facei];
                const scalarList& weights = allWeights[facei];

                forAll(slots, i)
                {
                    cop(result[facei], facei, fld[slots[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& slots = allSlots[facei];
            const scalarList& weights = allWeights[facei];

            forAll(slots, i)
            {
                cop(result[facei], facei, fld[slots[i]], weights[i]);
            }
        }
    }
*/
    foamExecutor exec;
    auto result_p = result.begin();
    //auto result_p = gpu_result.begin();
    const auto fld_p = fld.cbegin();
    const auto allSlots_p = allSlots.cbegin();
    const auto allWeights_p = allWeights.cbegin();
    auto defaultValues_p = defaultValues.begin();
    const auto weightsSum_p = weightsSum.cbegin();
    const auto flatSlots_p = flatSlots.cbegin();
    const auto flatWeights_p = flatWeights.cbegin();
    const auto flatIdx_p = flatIdx.cbegin();

/*
    Info<< "device valid fld: "<<isDeviceValid(fld_p)<<endl;
    Info<< "device valid defaultValues: "<<isDeviceValid(defaultValues_p)<<endl;
    Info<< "device valid result: "<<isDeviceValid(result_p)<<endl;
    Info<< "device valid flatSlots: "<<isDeviceValid(flatSlots_p)<<endl;
    Info<< "device valid flatWeights: "<<isDeviceValid(flatWeights_p)<<endl;
    Info<< "device valid flatIdx: "<<isDeviceValid(flatIdx_p)<<endl;
    Info<< "device valid allSlots: "<<isDeviceValid(allSlots_p)<<endl;
    Info<< "device valid allWeights: "<<isDeviceValid(allWeights_p)<<endl;
*/

    if (lowWeightCorrection > 0)
    {
        auto Lambda = [=](label facei){
            if (weightsSum_p[facei] < lowWeightCorrection)
            {
                result_p[facei] = defaultValues_p[facei];
            }
        };
        auto Lambda2 = [=](label facei){
            if (weightsSum_p[facei] >= lowWeightCorrection)
            {
                for (size_t i = flatIdx_p[facei]; i < flatIdx_p[facei+1]; i++)
                {
                    cop(result_p[facei], facei, fld_p[flatSlots_p[i]], flatWeights_p[i]);
                }
            }
        };

        exec.parallelFor(Lambda,result.size());
        exec.parallelFor(Lambda2,result.size());
    }
    else
    {
        auto Lambda2 = [=](label facei){
            for (size_t i = flatIdx_p[facei]; i < flatIdx_p[facei+1]; i++)
            {
                cop(result_p[facei], facei, fld_p[flatSlots_p[i]], flatWeights_p[i]);
            };
        };
        exec.parallelFor(Lambda2,result.size());
    }

/*
    // check result
    if constexpr(std::is_same<Type,vector>::value || std::is_same<Type,scalar>::value ){
        equal = true;
        forAll(result,facei)
        {
            if(mag(result[facei] - gpu_result[facei]) > 1e-6){
                equal = false;
                break;
            }
        }
        Info << "results equal: "<<equal<<endl;
        //Info <<"results:" <<  result <<endl;
        //Info <<"gpu results:" <<  gpu_result <<endl;
    }
*/
}


template<class Type>
void Foam::AMIInterpolation::weightedSum
(
    const bool interpolateToSource,
    const UList<Type>& fld,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    weightedSum
    (
        lowWeightCorrection_,
        (interpolateToSource ? srcAddress_ : tgtAddress_),
        (interpolateToSource ? srcWeights_ : tgtWeights_),
        (interpolateToSource ? flatSrcAddress_ : flatTgtAddress_),
        (interpolateToSource ? flatSrcWeights_ : flatTgtWeights_),
        (interpolateToSource ? flatSrcIdx_ : flatTgtIdx_),
        (interpolateToSource ? srcWeightsSum_ : tgtWeightsSum_),
        fld,
        multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
        result,
        defaultValues
    );
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToTarget");

    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }
    else if
    (
        (lowWeightCorrection_ > 0)
     && (defaultValues.size() != tgtAddress_.size())
    )
    {
        FatalErrorInFunction
            << "Employing default values when sum of weights falls below "
            << lowWeightCorrection_
            << " but supplied default field size is not equal to target "
            << "patch size" << nl
            << "    default values = " << defaultValues.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << abort(FatalError);
    }

    result.setSize(tgtAddress_.size());
    List<Type> work(poolSwitch(true));

    if (distributed())
    {
        const mapDistribute& map = srcMapPtr_();
        work.resize_nocopy(map.constructSize());
        SubList<Type>(work, fld.size()) = fld;  // deep copy
        map.distribute(work);
    }


    if constexpr(std::is_same<CombineOp,multiplyWeightedOp<Type, plusEqOp<Type>>>::value){
        weightedSum
        (
            lowWeightCorrection_,
            tgtAddress_,
            tgtWeights_,
            flatTgtAddress_,
            flatTgtWeights_,
            flatTgtIdx_,
            tgtWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }else{
        weightedSum
        (
            lowWeightCorrection_,
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToSource");

    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }
    else if
    (
        (lowWeightCorrection_ > 0)
     && (defaultValues.size() != srcAddress_.size())
    )
    {
        FatalErrorInFunction
            << "Employing default values when sum of weights falls below "
            << lowWeightCorrection_
            << " but number of default values is not equal to source "
            << "patch size" << nl
            << "    default values = " << defaultValues.size() << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << abort(FatalError);
    }

    result.setSize(srcAddress_.size());
    List<Type> work(poolSwitch(true));

    if (distributed())
    {
        const mapDistribute& map = tgtMapPtr_();
        work.resize_nocopy(map.constructSize());
        SubList<Type>(work, fld.size()) = fld;  // deep copy
        map.distribute(work);
    }

    if constexpr(std::is_same<CombineOp,multiplyWeightedOp<Type, plusEqOp<Type>>>::value){
        weightedSum
        (
            lowWeightCorrection_,
            srcAddress_,
            srcWeights_,
            flatSrcAddress_,
            flatSrcWeights_,
            flatSrcIdx_,
            srcWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }else{
        weightedSum
        (
            lowWeightCorrection_,
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(srcAddress_.size(), Zero);

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(tgtAddress_.size(), Zero);

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


// ************************************************************************* //
