#pragma once

namespace CESDSOL
{
	template<typename MatrixType, typename VectorType>
	class LinearSolver
	{
	public:
		virtual bool Solve(const MatrixType& matrix, const VectorType& y, VectorType& x) = 0;
		virtual ~LinearSolver() = default;
	};
}