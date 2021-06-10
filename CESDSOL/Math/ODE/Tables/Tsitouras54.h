#pragma once

#include <array>

namespace CESDSOL
{
	struct Tsitouras54
	{
		static constexpr size_t AccuracyOrder = 5;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 6;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.161,
			0.327,
			0.9,
			0.9800255409045097,
			1
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.161, 0, 0, 0, 0, 0},
			{-0.008480655492356989,
	0.335480655492357, 0, 0, 0,
	0}, {2.897153057105494, -6.359448489975075, 4.362295432869582, 0, 0,
	 0}, {5.325864828439257, -11.74888356406283,
	7.495539342889837, -0.09249506636175525, 0,
	0}, {5.86145544294642, -12.92096931784711,
	8.159367898576159, -0.071584973281401, -0.02826905039406838,
	0}, {0.09646076681806523, 0.01, 0.4798896504144996,
	1.379008574103742, -3.290069515436081, 2.324710524099774}
		  } };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 4 };
		static constexpr std::array<std::array<double, StepCount + 1>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				-0.00178001105222577714,
				-0.0008164344596567469,
				0.007880878010261995,
				-0.1447110071732629,
				0.5823571654525552,
				-0.45808210592918697,
				0.015151515151515152
			}
		} };

		static constexpr bool IsDenseOutputSupported = true;
		static constexpr size_t InterpolationOrder = 4;
		static constexpr size_t DenseOutputStepCount = 0;
		static constexpr std::array<std::array<double, InterpolationOrder + 1>, StepCount + 1 + DenseOutputStepCount> DenseOutputCoefficients
		{ {{0, 1., -2.763706197274826, 2.913255461821913, -1.053088497729022}, {0, 0,
	0.1317, -0.2234, 0.1017}, {0, 0, 3.930296236894752, -5.941033872131505,
	 2.490627285651253}, {0, 0, -12.41107716693368,
	30.33818863028232, -16.5481028892449}, {0, 0,
	37.50931341651104, -88.1789048947664,
	47.37952196281928}, {0, 0, -27.89652628919729,
	65.09189467479366, -34.87065786149661}, {0, 0, 1.5, -4, 2.5}}
		};
	};
}