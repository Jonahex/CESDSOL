#pragma once

#include <array>

namespace CESDSOL
{
	struct Midpoint21
	{
		static constexpr size_t AccuracyOrder = 2;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 2;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0,
			0.5
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.5}, {0, 1}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 1 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{1, -1}
		} };

		static constexpr bool IsDenseOutputSupported = false;
	};
}