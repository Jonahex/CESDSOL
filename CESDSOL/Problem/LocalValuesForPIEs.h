#pragma once

#include "Math/LinearAlgebra.h"
#include "Problem/ProblemType.h"

#include <array>

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldTypeArg = double>
	class LocalValuesForPIEs
	{
	public:
		using FieldType = FieldTypeArg;

		std::array<CoordinateType, Dimension> Point;
		Array<FieldType> PIEValues;

		LocalValuesForPIEs(size_t pieCount) noexcept
			: PIEValues(pieCount)
		{}
	};
}