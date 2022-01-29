#pragma once

#include "Math/LinearAlgebra.h"
#include "Math/MKL/ILUPreconditioner.h"

namespace CESDSOL::MKL
{
	class ILUT final
		: public ILUPreconditioner<CESDSOL::CSRMatrix<double>, Vector<double>>
	{
	public:
		ILUT() noexcept;

		bool Setup(const CESDSOL::CSRMatrix<double>& matrix, const Vector<double>& y) noexcept override;

		void SetNormalizeZeroDiagonal(bool value) noexcept;
		void SetZeroDiagonalNormalizer(double value) noexcept;

	private:
		std::array<MKL_INT, 128> intParameters;
		std::array<double, 128> doubleParameters;

		MakeProperty(maxfil, MaxFillIn, MKL_INT, 5);
		MakeProperty(tolerance, ThresholdTolerance, double, 1e-10);
	};
}