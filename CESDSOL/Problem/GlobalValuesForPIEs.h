#pragma once

#include "Math/LinearAlgebra.h"
#include "Problem/ProblemType.h"

namespace CESDSOL
{
	template<ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class GlobalValuesForPIEs
	{
	public:
		const Array<FieldType>& PIEValues;

		GlobalValuesForPIEs(const Array<FieldType>& pieValues) noexcept
			: PIEValues(pieValues)
		{}
	};
}