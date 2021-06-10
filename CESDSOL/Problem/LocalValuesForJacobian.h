#pragma once

#include "Problem/LocalValues.h"

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class LocalValuesForJacobian
		: public LocalValues<Dimension, Type, CoordinateType, FieldType>
	{
	public:
		ThreeLevelArray<FieldType> LVDEJacobianComponentValues;
		ThreeLevelArray<FieldType> ReductionJacobianComponentValues;

		LocalValuesForJacobian(
			size_t pieCount,
			size_t vieCount,
			size_t fieldsCount,
			LevelStructure<2> derivativesStructure,
			size_t variableDependentExpressionCount,
			LevelStructure<3> lvdeJacobianComponentValues,
			LevelStructure<3> reductionJacobianComponentValues
		) noexcept
			: LocalValues<Dimension, Type, CoordinateType, FieldType>(
				pieCount,
				vieCount,
				fieldsCount,
				derivativesStructure,
				variableDependentExpressionCount)
			, LVDEJacobianComponentValues(lvdeJacobianComponentValues)
			, ReductionJacobianComponentValues(reductionJacobianComponentValues)
		{}
	};
}