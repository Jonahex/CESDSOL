#pragma once

#include "Problem/GlobalValues.h"

namespace CESDSOL
{
	template<ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class GlobalValuesForJacobian
		: public GlobalValues<Type, CoordinateType, FieldType>
	{
	public:
		const TwoLevelArray<FieldType>& GVDEJacobianComponentValues;

		GlobalValuesForJacobian(
			const Array<FieldType>& pieValues,
			const Array<FieldType>& parameters,
			const Array<FieldType>& vieValues,
			const std::span<FieldType>& discreteVariables,
			const Array<FieldType>& vdeValues,
			const Array<FieldType>& reductionValues,
			const TwoLevelArray<FieldType>& gvdeJacobianComponentValues
		) noexcept
			: GlobalValues<Type, CoordinateType, FieldType>
			( pieValues
			, parameters
			, vieValues
			, discreteVariables
			, vdeValues
			, reductionValues
			)
			, GVDEJacobianComponentValues(gvdeJacobianComponentValues)
		{}
	};
}