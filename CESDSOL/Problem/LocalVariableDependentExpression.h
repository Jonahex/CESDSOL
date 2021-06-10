#pragma once

#include "Problem/GlobalValuesForJacobian.h"
#include "Problem/LocalValuesForJacobian.h"

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class LocalVariableDependentExpression
	{
	public:
		using CurrentLocalValues = LocalValues<Dimension, Type, CoordinateType, FieldType>;
		using CurrentLocalValuesForJacobian = LocalValuesForJacobian<Dimension, Type, CoordinateType, FieldType>;
		using CurrentGlobalValues = GlobalValues<Type, CoordinateType, FieldType>;
		using CurrentGlobalValuesForJacobian = GlobalValuesForJacobian<Type, CoordinateType, FieldType>;

		std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> Expression;
		TwoLevelArray<std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>> JacobianComponents;

		[[nodiscard]] bool IsValid() const noexcept
		{
			return static_cast<bool>(Expression);
		}
	};
}