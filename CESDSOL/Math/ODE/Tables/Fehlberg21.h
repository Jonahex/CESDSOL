#pragma once

#include <array>

namespace CESDSOL
{
	struct Fehlberg21
	{
		static constexpr size_t AccuracyOrder = 2;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 3;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0,
			0.5,
			1.
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.5},
			{0.00390625, 0.99609375},
			{0.001953125, 0.99609375, 0.001953125}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 1 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{0.00390625, 0., -0.00390625}
		} };

		static constexpr bool IsDenseOutputSupported = false;
	};
}