#pragma once

#include <array>

namespace CESDSOL
{
	struct BogackiShampine32
	{
		static constexpr size_t AccuracyOrder = 3;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 4;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.5,
			0.75,
			1
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{ 0.5 },
			{ 0., 0.75 },
			{0.2222222222222222, 0.3333333333333333, 0.4444444444444444},
			{0.2222222222222222, 0.3333333333333333, 0.4444444444444444}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 2 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				0.06944444444444445,
				-0.08333333333333333,
				-0.1111111111111111,
				0.125
			}
		} };

		static constexpr bool IsDenseOutputSupported = false;
	};
}
