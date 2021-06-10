#pragma once

#include "Problem/LocalValuesForPIEs.h"

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class LocalValuesForVIEs
		: public LocalValuesForPIEs<Dimension, Type, CoordinateType, FieldType>
	{
	public:
		Array<FieldType> VIEValues;

		LocalValuesForVIEs(
			size_t pieCount,
			size_t vieCount
		) noexcept
			: LocalValuesForPIEs<Dimension, Type, CoordinateType, FieldType>(pieCount)
			, VIEValues(vieCount)
		{}
	};
}