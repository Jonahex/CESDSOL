#pragma once

#include "Math/NonlinearSolver.h"
#include "Utils/Aliases.h"

namespace CESDSOL
{
	template<typename ProblemType>
	class ParametricSweeper
	{
	public:
		using ValueType = typename ProblemType::FieldType;

		struct OutputInfo
		{
			bool success = true;
		};

		virtual uptr<OutputInfo> Sweep() const noexcept = 0;
		virtual ~ParametricSweeper() = default;
	
	protected:
		ParametricSweeper
			( sptr<ProblemType> aProblem
			, uptr<NonlinearSolver<ProblemType>> aSolver
			) noexcept
			: problem(std::move(aProblem))
			, solver(std::move(aSolver))
		{}

		sptr<ProblemType> problem;
		uptr<NonlinearSolver<ProblemType>> solver;
	};
}