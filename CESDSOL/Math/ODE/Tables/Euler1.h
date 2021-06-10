#pragma once

#include <array>

namespace CESDSOL
{
	struct Euler1
	{
		static constexpr size_t AccuracyOrder = 1;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 1;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{1}
		} };

		static constexpr bool IsAdaptive = false;

		static constexpr bool IsDenseOutputSupported = false;
	};
}