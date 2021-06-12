#pragma once

#include "Math/LinearAlgebra.h"
#include "Math/Preconditioner.h"

namespace CESDSOL::MKL
{
	class ILUT final
		: public Preconditioner<CESDSOL::CSRMatrix<double>>
	{
	public:
		ILUT() noexcept;

		uptr<OutputInfo> Apply(CSRMatrix<double>& matrix) const noexcept override;

		void SetNormalizeZeroDiagonal(bool value) noexcept;
		void SetZeroDiagonalNormalizer(double value) noexcept;

	private:
		std::array<MKL_INT, 128> intParameters;
		std::array<double, 128> doubleParameters;

		MakeProperty(maxfil, MaxFillIn, MKL_INT, 5);
		MakeProperty(tolerance, ThresholdTolerance, double, 1e-10);
	};
}