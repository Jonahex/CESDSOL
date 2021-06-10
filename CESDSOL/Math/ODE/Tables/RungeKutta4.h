#pragma once

#include <array>

namespace CESDSOL
{
	struct RungeKutta4
	{
		static constexpr size_t AccuracyOrder = 4;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 4;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.5,
			0.5,
			1.
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{ 0.5 },
			{ 0., 0.5 },
			{ 0., 0., 1. },
			{ 1. / 6., 1. / 3., 1. / 3., 1. / 6. }
		} };

		static constexpr bool IsAdaptive = false;

		static constexpr bool IsDenseOutputSupported = false;
	};
}