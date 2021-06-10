#pragma once

#include "Problem/GlobalValuesForJacobian.h"

namespace CESDSOL
{
	template<ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class GlobalVariableDependentExpression
	{
	public:
		using CurrentGlobalValues = GlobalValues<Type, CoordinateType, FieldType>;
		using CurrentGlobalValuesForJacobian = GlobalValuesForJacobian<Type, CoordinateType, FieldType>;

		std::function<FieldType(const CurrentGlobalValues&)> Expression;
		Array<std::function<FieldType(const CurrentGlobalValuesForJacobian&)>> JacobianComponents;

		[[nodiscard]] bool IsValid() const noexcept
		{
			return static_cast<bool>(Expression);
		}
	};
}