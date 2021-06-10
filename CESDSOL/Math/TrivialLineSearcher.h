#pragma once

#include "Math/LineSearcher.h"
#include "Utils/Utils.h"

namespace CESDSOL
{
	template<typename ProblemType>
	class TrivialLineSearcher final
		: public LineSearcher<ProblemType>
	{
	public:
		using LineSearcher<ProblemType>::ValueType;
		using LineSearcher<ProblemType>::OutputInfo;
		
		MakeProperty(shiftFactor, ShiftFactor, double, 1)

		TrivialLineSearcher(double aShiftFactor = 1) noexcept
			: shiftFactor(aShiftFactor)
		{}

		uptr<OutputInfo> Solve(ProblemType& problem, const Vector<ValueType>& shift) const noexcept override
		{
			AXPY(shiftFactor, shift, problem.GetVariables().Flatten());
			problem.SetVariablesUpdated();
			return std::make_unique<OutputInfo>();
		}
	};
}