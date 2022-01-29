#pragma once

#include "Math/LinearSolver.h"

namespace CESDSOL
{
	template<typename MatrixType, typename VectorType>
	class Preconditioner : public LinearSolver<MatrixType, VectorType>
	{
	public:
		virtual bool Setup(const MatrixType& matrix, const VectorType& y) noexcept 
		{
			return true;
		};
	};
}