#pragma once

#include "Problem/GlobalValuesForPIEs.h"

namespace CESDSOL
{
	template<ProblemType Type = ProblemType::Stationary, typename CoordinateType = double, typename FieldType = double>
	class GlobalValuesForVIEs
		: public GlobalValuesForPIEs<Type, CoordinateType, FieldType>
	{
	public:
		const Array<FieldType>& Parameters;
		const Array<FieldType>& VIEValues;

		GlobalValuesForVIEs(
			const Array<FieldType>& pieValues,
			const Array<FieldType>& parameters,
			const Array<FieldType>& vieValues
		) noexcept
			: GlobalValuesForPIEs<Type, CoordinateType, FieldType>(pieValues)
			, Parameters(parameters)
			, VIEValues(vieValues)
		{}
	};

	template<ProblemType Type, typename CoordinateType, typename FieldType>
	requires (Type == ProblemType::ExplicitTransient || Type == ProblemType::ImplicitTransient)
	class GlobalValuesForVIEs<Type, CoordinateType, FieldType>
		: public GlobalValuesForVIEs<ProblemType::Stationary, CoordinateType, FieldType>
	{
	private:
		using BaseType = GlobalValuesForVIEs<ProblemType::Stationary, CoordinateType, FieldType>;
	
	public:
		CoordinateType Time;

		using BaseType::BaseType;
	};
}