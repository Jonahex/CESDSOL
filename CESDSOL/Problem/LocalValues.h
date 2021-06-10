#pragma once

#include "Problem/LocalValuesForVIEs.h"

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class LocalValues
		: public LocalValuesForVIEs<Dimension, Type, CoordinateType, FieldType>
	{
	public:
		Array<FieldType> FieldValues;
		TwoLevelArray<FieldType> DerivativeValues;
		Array<FieldType> VDEValues;
		FieldType IntegrationWeight;

		LocalValues(
			size_t pieCount,
			size_t vieCount,
			size_t fieldsCount,
			LevelStructure<2> derivativesStructure,
			size_t variableDependentExpressionCount
		) noexcept
			: LocalValuesForVIEs<Dimension, Type, CoordinateType, FieldType>(pieCount, vieCount)
			, FieldValues(fieldsCount)
			, DerivativeValues(derivativesStructure)
			, VDEValues(variableDependentExpressionCount)
		{}
	};
}