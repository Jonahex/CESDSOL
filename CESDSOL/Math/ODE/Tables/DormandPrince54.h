#pragma once

#include <array>

namespace CESDSOL
{
	struct DormandPrince54
	{
		static constexpr size_t AccuracyOrder = 5;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 6;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.2,
			0.3,
			0.8,
			0.8888888888888888,
			1
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{},
			{0.2},
			{0.075, 0.225}, {0.9777777777777777, -3.733333333333333, 3.555555555555555},
			{2.952598689224204, -11.59579332418839, 9.822892851699436, -0.2908093278463649},
			{2.846275252525253, -10.75757575757576, 8.906422717743473, 0.2784090909090909, -0.2735313036020583},
			{0.09114583333333333, 0, 0.4492362982929021, 0.6510416666666666, -0.322376179245283, 0.130952380952381}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 4 };
		static constexpr std::array<std::array<double, StepCount + 1>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				-0.0012326388888888888,
				0,
				0.0042527702905061394,
				-0.03697916666666667,
				0.05086379716981132,
				-0.0419047619047619,
				0.025
			}
		} };

		static constexpr bool IsDenseOutputSupported = true;
		static constexpr size_t DenseOutputStepCount = 0;
		static constexpr size_t InterpolationOrder = 4;
		static constexpr std::array<std::array<double, InterpolationOrder + 1>, StepCount + 1 + DenseOutputStepCount> DenseOutputCoefficients =
		{{
			{0, 1., -2.853580065386284, 3.0717434641059, -1.127017565386284},
			{},
			{0, 0, 4.023133379230305, -6.249321565289, 2.675424484351598},
			{0, 0, -3.732401961588505, 10.06897058984368, -5.685526961588504},
			{0, 0, 2.554803830184942, -6.399112377351017, 3.521932367920791},
			{0, 0, -1.374424114218603, 3.272657752246729, -1.767281257075746},
			{0, 0, 1.382468931778144, -3.764937863556288, 2.382468931778144}
		} };
	};
}
