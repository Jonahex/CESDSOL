#pragma once

#include <array>

namespace CESDSOL
{
	struct OwrenZennaro32
	{
		static constexpr size_t AccuracyOrder = 3;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 3;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.5217391304347826,
			0.8
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{ 0.5217391304347826 },
			{ -0.18133333333333335, 0.9813333333333333 },
			{0.2152777777777778, 0.4592013888888889, 0.3255208333333333}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 2 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				-0.1736111111111111,
				0.4991319444444444,
				-0.3255208333333333
			}
		} };

		static constexpr bool IsDenseOutputSupported = true;
		static constexpr size_t InterpolationOrder = 3;
		static constexpr size_t DenseOutputStepCount = 0;
		static constexpr std::array<std::array<double, InterpolationOrder + 1>, StepCount + 1 + DenseOutputStepCount> DenseOutputCoefficients
		{ {
			{0., 1., -1.3541666666666667, 0.5694444444444444},
			{0., 0., 1.3776041666666667, -0.9184027777777778},
			{0., 0., 0.9765625, -0.6510416666666666},
			{0., 0., -1., 1.}
		} };
	};
}