#pragma once

#include <array>

namespace CESDSOL
{
	struct Ralston21
	{
		static constexpr size_t AccuracyOrder = 2;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 2;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0,
			0.6666666666666666
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.6666666666666666}, {0.25, 0.75}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 1 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{0.6666666666666666, -0.6666666666666666}
		} };

		static constexpr bool IsDenseOutputSupported = false;
	};
}