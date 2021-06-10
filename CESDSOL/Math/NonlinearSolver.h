#pragma once

#include "Utils/Aliases.h"

namespace CESDSOL
{
	template<typename ProblemType>
	class NonlinearSolver
	{
	public:
		struct OutputInfo
		{
			bool success = true;
		};

		virtual uptr<OutputInfo> Solve(ProblemType& problem) const noexcept = 0;
		virtual ~NonlinearSolver() = default;
	};

	template<template<typename> typename SolverTemplate, typename ProblemType, typename... CtorArgTypes>
	auto MakeNonlinearSolver(const ProblemType& problem, CtorArgTypes&&... args)
	{
		return std::make_unique<SolverTemplate<ProblemType>>(std::forward<CtorArgTypes>(args)...);
	}
}