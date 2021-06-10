#pragma once

#include <array>

namespace CESDSOL
{
	struct OwrenZennaro54
	{
		static constexpr size_t AccuracyOrder = 5;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 7;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.16666666666666666,
			0.25,
			0.5,
			0.5,
			0.6428571428571429,
			0.875
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.1666666666666667 },
			{0.0625, 0.1875 },
			{0.25, -0.75 },
			{-0.75, 3.75, -3., 0.5 },
			{0.2689504373177842, -0.7084548104956269, 0.8658892128279884,
  0.1546230737192836, 0.06184922948771345 },
			{-0.02947695035460993, 0.1850066489361702,
  0.4802345261121857, -0.5337849069148937, -0.01309009308510638,
  0.7861107753062541 },
			{0.08783068783068783, 0, 0.3006060606060606,
   0.2277777777777778, 0.02777777777777778, 0.06218596218596219,
  0.2938217338217338}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 4 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				-0.19894179894179895,
				0.,
				0.9115151515151515,
				-1.9777777777777779,
				-0.1111111111111111,
				1.67013727013727,
				-0.2938217338217338
			}
		} };

		static constexpr bool IsDenseOutputSupported = true;
		static constexpr size_t InterpolationOrder = 5;
		static constexpr size_t DenseOutputStepCount = 0;
		static constexpr std::array<std::array<double, InterpolationOrder + 1>, StepCount + 1 + DenseOutputStepCount> DenseOutputCoefficients
		{ {{0, 1, -4.01953601953602, 7.282458282458283, -6.067155067155067,
  1.892063492063492}, {}, {0, 0,
  7.14965034965035, -20.31142191142191,
  20.67692307692308, -7.214545454545455}, {0, 0, -2.365384615384615,
  13.50854700854701, -18.78205128205128,
  7.866666666666666}, {0, 0, -1.211538461538462,
  4.534188034188034, -5.294871794871795, 2}, {0, 0, -1.219801565955412,
  1.195883888191581, 1.578566732412886, -1.492463092463092}, {0, 0,
  3.051225697379544, -11.97888607119376,
  16.27320371935756, -7.051721611721612}, {0, 0, -1.384615384615385,
  5.769230769230769, -8.384615384615385, 4}} };
	};
}