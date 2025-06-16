/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

    Copyright (C) 2022 Stefano Oliani
    Copyright (C) 2025 Cineca
-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM & OpenFOAM.

    ICSFOAM/OpenFOAM are free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM/OpenFOAM are distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM/OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::expandData(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

	const scalar myAngle = 360.0/scalar(ncp);

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField.ref();

    label pfSize = pf.size();

    foamExecutor exec;
    auto expandFieldPtr = expandField.begin();
    const auto pfPtr = pf.cbegin();
    const auto curRotationsPtr = curRotations_.cbegin();

	auto Lambda = [=](label faceI)
	{
	    for (label copyI = 0; copyI < ncp; copyI++)
	    {
	         const label offset = copyI*pfSize;
		     const label zId = this->whichFace(this->start() + faceI);
		     const tensor& curRotation = curRotationsPtr[copyI];
		     expandFieldPtr[offset + zId] = Foam::transform(curRotation, pfPtr[faceI]);
	    }
	};
	exec.parallelFor(Lambda,pf.size());

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::expandData(const Field<Type>& pf, label cmpt) const
{
	Info<<"should not call it"<<endl;

    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

	//const scalar myAngle = 360.0/scalar(ncp);

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    //Field<Type>& expandField = texpandField.ref();

    //foamExecutor exec;
    //auto expandFieldPtr = expandField.begin();
    //const auto pfPtr = pf.cbegin();

    //for (label copyI = 0; copyI < ncp; copyI++)
    //{
    //	// Calculate transform
	//	const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

	//	const label offset = copyI*pf.size();

	//	auto Lambda = [=](label faceI)
	//	{
	//		 const label zId = this->whichFace(this->start() + faceI);
	//		 expandFieldPtr[offset + zId] = Foam::transform(curRotation, pfPtr[faceI]);
	//	};

	//	exec.parallelFor(Lambda,pf.size());
	//	//forAll (pf, faceI)
	//	//{
	//	//	 const label zId = this->whichFace(this->start() + faceI);
	//	//	 expandField[offset + zId] = Foam::transform(curRotation, pf[faceI]);
	//	//}
    //}

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::untransfExpandData(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField.ref();

    label pfSize = pf.size();

    foamExecutor exec;
    auto expandFieldPtr = expandField.begin();
    const auto pfPtr = pf.cbegin();

    auto Lambda = [=](label faceI)
	{
        for (label copyI = 0; copyI < ncp; copyI++)
        {
	    	const label offset = copyI*pfSize;
		    const label zId = this->whichFace(this->start() + faceI);
		    expandFieldPtr[offset + zId] = pfPtr[faceI];
        }
	};

	exec.parallelFor(Lambda,pf.size());

    return texpandField;
}


template<class Type>
Foam::UList<Type>
Foam::overlapAMIPolyPatch::expandData(const UList<Type>& defaultValues) const
{
	if (defaultValues.size())
	{
		// Check and expand the field from patch size to zone size
		if (defaultValues.size() != this->size())
		{
			FatalErrorIn
			(
				"tmp<Field<Type> > overlapAMIPolyPatch::expandData"
				"("
				"    UList<Type>& defaultValues"
				") const"
			)   << "Incorrect patch field size.  Field size: "
				<< defaultValues.size() << " patch size: " << this->size()
				<< abort(FatalError);
		}

		const label ncp = nCopies();

		const scalar myAngle = 360.0/scalar(ncp);

		Type dfl = *defaultValues.cdata();
		Type* dflPtr = &dfl;

		UList<Type> expandField(dflPtr, nCopies_*defaultValues.size(),true);

        label dfSize = defaultValues.size();

        foamExecutor exec;
        auto expandFieldPtr = expandField.begin();
        const auto defaultValuesPtr = defaultValues.cbegin();
        const auto curRotationsPtr = curRotations_.cbegin();

		auto Lambda = [=](label faceI)
		{
		    for (label copyI = 0; copyI < ncp; copyI++)
		    {
		        const label offset = copyI*dfSize;
			     const label zId = this->whichFace(this->start() + faceI);
			     const tensor& curRotation = curRotationsPtr[copyI];
			     expandFieldPtr[offset + zId] = Foam::transform(curRotation, defaultValuesPtr[faceI]);
		    }
		};

		exec.parallelFor(Lambda,defaultValues.size());

		return expandField;
	}

	return defaultValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    nvtxRangePushA("oAMIfvpoly_int1");
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type>  expandDefault = neighbPatch().expandData(defaultValues);

	tmp<Field<Type> > tresult(new Field<Type>());
	Field<Type>& result = tresult.ref();
    nvtxRangePop();

    if (owner())
    {
    nvtxRangePushA("oAMIfvpoly_int1A1");
        result = AMI().interpolateToSource(expandData, expandDefault);
nvtxRangePop();
        // Truncate to size
    nvtxRangePushA("oAMIfvpoly_int1A2");
        result.setSize(this->size());

nvtxRangePop();
        return tresult;
    }
    else
    {
    nvtxRangePushA("oAMIfvpoly_int1B");
        result = neighbPatch().AMI().interpolateToTarget(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

nvtxRangePop();
        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    nvtxRangePushA("oAMIfvpoly_int2");
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld);

	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type> expandDefault = neighbPatch().expandData(defaultValues);

    if (owner())
    {
        AMI().interpolateToSource
        (
            expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
        	expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
nvtxRangePop();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
	label cmpt,
    const UList<Type>& defaultValues
) const
{
    nvtxRangePushA("oAMIfvpoly_int3");
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, cmpt);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type>  expandDefault = neighbPatch().expandData(defaultValues);

	tmp<Field<Type> > tresult(new Field<Type>());
	Field<Type>& result = tresult.ref();

    if (owner())
    {
        result = AMI().interpolateToSource(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

nvtxRangePop();
        return tresult;
    }
    else
    {
        result = neighbPatch().AMI().interpolateToTarget(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

nvtxRangePop();
        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
	label cmpt,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), cmpt, defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
	label cmpt,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    nvtxRangePushA("oAMIfvpoly_int4");
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, cmpt);

	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type> expandDefault = neighbPatch().expandData(defaultValues);

    if (owner())
    {
        AMI().interpolateToSource
        (
            expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
        	expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
nvtxRangePop();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::untransfInterp
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().untransfExpandData(fld);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type>  expandDefault = neighbPatch().expandData(defaultValues);

	tmp<Field<Type> > tresult(new Field<Type>());
	Field<Type>& result = tresult.ref();

    if (owner())
    {
        result = AMI().interpolateToSource(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
    else
    {
        result = neighbPatch().AMI().interpolateToTarget(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::untransfInterp
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return untransfInterp(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::untransfInterp
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().untransfExpandData(fld);

	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type> expandDefault = neighbPatch().expandData(defaultValues);

    if (owner())
    {
        AMI().interpolateToSource
        (
            expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
        	expandData,
            cop,
            result,
			expandDefault
        );

        // Truncate to size
		result.setSize(this->size());
    }
}

// ************************************************************************* //
