#pragma once

#include <array>

namespace CESDSOL
{
	struct OwrenZennaro43
	{
		static constexpr size_t AccuracyOrder = 4;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 5;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0.,
			0.16666666666666666,
			0.2972972972972973,
			0.6470588235294118,
			0.8666666666666667
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{{
			{},
			{0.1666666666666667},
			{0.03214024835646457, 0.26515704894083280},
			{0.6895990230002036, -1.699369020964787, 1.656828821493996},
			{-0.09002509947964493, 0.6817777777777778, -0.2402791551882461, 0.5151931435567799},
			{0.08990252172070354, 0, 0.4360623278236915, 0.1842858372687918, 0.2897493131868132}
		}};

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 3 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				0.18833439287984743,
				0.,
				-0.5303460743801653,
				0.6317609946871311,
				-0.2897493131868132
			}
		} };

		static constexpr bool IsDenseOutputSupported = true;
		static constexpr size_t InterpolationOrder = 4;
		static constexpr size_t DenseOutputStepCount = 0;
		static constexpr std::array<std::array<double, InterpolationOrder + 1>, StepCount + 1 + DenseOutputStepCount> DenseOutputCoefficients
		{ {
			{0, 1, -2.781642022100038, 2.922894131082889, -1.051349587262148},
				{},
				{0, 0, 3.734823907009022, -5.725398502723277, 2.426636923537947},
			{0, 0, -0.2176560796074155, 1.172455508289998, -0.7705135914137909},
			{0, 0, -1.995067790034393, 5.149132832816039, -2.864315729594833},
			{0, 0, 1.259541984732824, -3.519083969465649, 2.259541984732825}
		} };
	};
}