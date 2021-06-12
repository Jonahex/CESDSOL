#pragma once

#include "Math/LinearAlgebra.h"
#include "Math/Preconditioner.h"

namespace CESDSOL::MKL
{
	class ILU0 final
		: public Preconditioner<CESDSOL::CSRMatrix<double>>
	{
	public:
		ILU0() noexcept;
		
		uptr<OutputInfo> Apply(CSRMatrix<double>& matrix) const noexcept override;

		void SetNormalizeZeroDiagonal(bool value) noexcept;
		void SetZeroDiagonalThreshold(double value) noexcept;
		void SetZeroDiagonalNormalizer(double value) noexcept;

	private:
		std::array<MKL_INT, 128> intParameters;
		std::array<double, 128> doubleParameters;
	};
}