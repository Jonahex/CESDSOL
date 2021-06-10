#pragma once

#include "Grid/Grid.h"
#include "Problem/BaseProblemDescriptor.h"

namespace CESDSOL
{
	template<size_t DimensionArg, template<typename> typename MatrixTypeArg = CSRMatrix, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class ExplicitTransientProblem;

	template<size_t DimensionArg, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class ExplicitTransientProblemDescriptor
		: public BaseProblemDescriptor<ProblemType::ExplicitTransient, DimensionArg, CoordinateTypeArg, FieldTypeArg>
	{
	private:
		using BaseType = BaseProblemDescriptor<ProblemType::ExplicitTransient, DimensionArg, CoordinateTypeArg, FieldTypeArg>;

	public:
		using typename BaseType::CoordinateType;
		using typename BaseType::FieldType;
		using typename BaseType::GridDescriptorType;
		using typename BaseType::CurrentLocalValues;
		using typename BaseType::CurrentGlobalValues;

		using BaseType::BaseType;

		static constexpr size_t Dimension = BaseType::Dimension;
		static constexpr ProblemType ProblemType = BaseType::ProblemType;

		template<template<typename> typename MatrixType = CSRMatrix>
		[[nodiscard]] auto
			MakeProblem(sptr<Grid<Dimension, CoordinateType>> grid, uptr<Discretization<Dimension, MatrixType, CoordinateType>> discretizer) const noexcept;
	};
}