#pragma once

#include "Grid/Grid.h"
#include "Math/LinearAlgebra.h"

namespace CESDSOL::ProblemUtils
{
	template<typename FieldType, typename GridType, typename ProblemType>
	void RebaseFields(const Vector<FieldType>& variables, const GridType& originalGrid, ProblemType& problem) noexcept
	{
		const auto fieldCount = problem.GetDescriptor().ContinuousEquationCount();
		const auto& newGrid = problem.GetGrid();
		const auto originalGridSize = originalGrid.GetSize();
		for (size_t gridIndex = 0; gridIndex < newGrid.GetSize(); ++gridIndex)
		{
			const auto weights = problem.GetDiscretizer().GetInterpolationWeightsVector(originalGrid, newGrid.GetCoordinates(gridIndex));
			for (size_t fieldIndex = 0; fieldIndex < fieldCount; ++fieldIndex)
			{
				const auto field = std::span(variables.begin() + fieldIndex * originalGridSize, originalGridSize);
				problem.SetVariable(fieldIndex, gridIndex, DotProduct(field, weights));
			}
		}
		for (size_t varIndex = 0; varIndex < problem.GetDescriptor().DiscreteEquationCount(); ++varIndex)
		{
			problem.SetVariable(fieldCount + varIndex, 0, variables[fieldCount * originalGridSize + varIndex]);
		}
	}

	template<typename ProblemType>
	void RebaseFields(const ProblemType& oldProblem, ProblemType& newProblem) noexcept
	{
		RebaseFields(oldProblem.GetVariables(), oldProblem.GetGrid(), newProblem);
	}
}