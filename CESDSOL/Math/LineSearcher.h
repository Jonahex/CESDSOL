#pragma once

#include "Math/LinearAlgebra.h"
#include "Utils/Aliases.h"

namespace CESDSOL
{
	template<typename ProblemType>
	class LineSearcher
	{
	public:
		using ValueType = typename ProblemType::FieldType;
		
		struct OutputInfo
		{
			bool success = true;
		};

		virtual uptr<OutputInfo> Solve(ProblemType& problem, const Vector<ValueType>& shift) const noexcept = 0;
		virtual ~LineSearcher() = default;
	};

	template<template<typename> typename SearcherTemplate, typename ProblemType, typename... CtorArgTypes>
	auto MakeLineSearcher(const ProblemType& problem, CtorArgTypes&&... args)
	{
		return std::make_unique<SearcherTemplate<ProblemType>>(std::forward<CtorArgTypes>(args)...);
	}
}