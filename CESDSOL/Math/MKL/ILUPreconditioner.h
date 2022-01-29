#pragma once

#include "Math/Preconditioner.h"

namespace CESDSOL::MKL
{
	template<typename MatrixType, typename VectorType>
	class ILUPreconditioner : public Preconditioner<MatrixType, VectorType>
	{
	public:
		bool Solve(const MatrixType& matrix, const VectorType& y, VectorType& x) override
		{
			TriangularSolve(preconditionedMatrix, y, temporaryArray, 1., false);
			TriangularSolve(preconditionedMatrix, temporaryArray, x, 1., true);
			return true;
		}

		bool Setup(const MatrixType& matrix, const VectorType& y) noexcept
		{
			temporaryArray = VectorType(y.size());
			return true;
		}

	protected:
		MatrixType preconditionedMatrix;

	private:
		VectorType temporaryArray;
	};
}