#pragma once

#include "Problem/GlobalValuesForJacobian.h"
#include "Problem/LocalValuesForJacobian.h"

namespace CESDSOL
{
	template<size_t Dimension, ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class Reduction
	{
	public:
		[[nodiscard]] static constexpr FieldType DefaultExternalFunction(FieldType x) noexcept
		{
			return x;
		}

		[[nodiscard]] static constexpr FieldType DefaultExternalJacobian(FieldType x) noexcept
		{
			return 1;
		}

		using CurrentLocalValues = LocalValues<Dimension, Type, CoordinateType, FieldType>;
		using CurrentGlobalValues = GlobalValues<Type, CoordinateType, FieldType>;
		using CurrentLocalValuesForJacobian = LocalValuesForJacobian<Dimension, Type, CoordinateType, FieldType>;
		using CurrentGlobalValuesForJacobian = GlobalValuesForJacobian<Type, CoordinateType, FieldType>;

		std::function<FieldType(FieldType)> ExternalFunction = DefaultExternalFunction;
		std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> InternalFunction;
		std::function<FieldType(FieldType)> ExternalFunctionJacobian = DefaultExternalJacobian;
		TwoLevelArray<std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>> InternalFunctionJacobianComponents;

		[[nodiscard]] bool IsValid() const noexcept
		{
			return ExternalFunction && InternalFunction;
		}
	};
}