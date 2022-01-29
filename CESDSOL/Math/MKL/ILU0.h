#pragma once

#include "Math/LinearAlgebra.h"
#include "Math/MKL/ILUPreconditioner.h"

namespace CESDSOL::MKL
{
	class ILU0 final
		: public ILUPreconditioner<CESDSOL::CSRMatrix<double>, Vector<double>>
	{
	public:
		ILU0() noexcept;
		
		bool Setup(const CESDSOL::CSRMatrix<double>& matrix, const Vector<double>& y) noexcept override;

		void SetNormalizeZeroDiagonal(bool value) noexcept;
		void SetZeroDiagonalThreshold(double value) noexcept;
		void SetZeroDiagonalNormalizer(double value) noexcept;

	private:
		std::array<MKL_INT, 128> intParameters;
		std::array<double, 128> doubleParameters;
	};
}