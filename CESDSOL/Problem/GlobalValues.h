#pragma once

#include "Problem/GlobalValuesForVIEs.h"

namespace CESDSOL
{
	template<ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class GlobalValues
		: public GlobalValuesForVIEs<Type, CoordinateType, FieldType>
	{
	public:
		std::span<FieldType> DiscreteVariables;
		const Array<FieldType>& VDEValues;
		const Array<FieldType>& ReductionValues;

		GlobalValues(
			const Array<FieldType>& pieValues,
			const Array<FieldType>& parameters,
			const Array<FieldType>& vieValues,
			const std::span<FieldType>& discreteVariables,
			const Array<FieldType>& vdeValues,
			const Array<FieldType>& reductionValues
		) noexcept
			: GlobalValuesForVIEs<Type, CoordinateType, FieldType>(pieValues, parameters, vieValues)
			, DiscreteVariables(discreteVariables)
			, VDEValues(vdeValues)
			, ReductionValues(reductionValues)
		{}
	};
}