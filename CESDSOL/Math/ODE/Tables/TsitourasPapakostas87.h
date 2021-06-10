#pragma once

#include <array>

namespace CESDSOL
{
	struct TsitourasPapakostas87
	{
		static constexpr size_t AccuracyOrder = 8;
		static constexpr bool IsExplicit = true;

		static constexpr size_t StepCount = 13;
		static constexpr std::array<double, StepCount> ButcherTableauFirstColumn =
		{
			0,
			0.06338028169014084,
  0.1027879458763643,
  0.15418191881454646,
  0.3875968992248062,
  0.4657534246575342,
  0.1554054054054054,
  1.0070921985815602,
  0.876141078561489,
  0.9120879120879121,
  0.959731543624161,
			1,
			1
		};
		static constexpr std::array<std::array<double, StepCount>, StepCount + 1> ButcherTableauMainPart =
		{ {
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.06338028169014084, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0}, {0.0194389804273365, 0.08334896544902781,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.03854547970363662, 0,
  0.1156364391109098, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.394365577701125,
  0, -1.481871932167337, 1.475103253691019, 0, 0, 0, 0, 0, 0, 0,
  0}, {0.0459944891076982, 0, 0, 0.2323507062639547,
  0.1874082292858813, 0, 0, 0, 0, 0, 0, 0}, {0.06005228953244051, 0,
  0, 0.1122038319463678, -0.03357232951906142, 0.01672161344565858, 0,
   0, 0, 0, 0, 0}, {-1.573329273208686, 0,
  0, -1.316708773022366, -11.72351529618177, 9.107825028173872,
  6.512820512820513, 0, 0, 0, 0, 0}, {-0.4810762562439125, 0,
  0, -6.65061036074639, -4.530206099782572, 3.894414525020157,
  8.634217645525526, 0.009401624788681498, 0, 0, 0,
  0}, {-0.7754121446230569, 0,
  0, -7.996604718235832, -6.726558607230182, 5.532184454327406,
  10.89757332024991, 0.0200916502800454, -0.03918604268037686, 0, 0,
  0}, {-1.189636324544999, 0,
  0, -7.128368483301214, -9.53722789710108, 7.574470108980868,
  11.26748638207092, 0.05100980122305832,
  0.08019413469508256, -0.1581961783984735, 0,
  0}, {-0.3920003904712727, 0, 0,
  3.916659042493856, -2.801745928908056,
  2.441204566481742, -2.418365577882472, -0.3394332629003293,
  0.1949645038310336, -0.1943717676250815, 0.5930888149805791,
  0}, {-1.484706308129189, 0,
  0, -2.390723588981498, -11.18430677284053, 8.720804667459817,
  7.33673830753461, 0.01289874999394761,
  0.0425832898426577, -0.05328834487981156, 0, 0},
			{0.04441161093250152, 0, 0, 0, 0,
  0.35395063113733116,
  0.2485219684184965,
  -0.3326913171720666,
  1.921248828652836,
  -2.7317783000882523,
  1.4012004409899175,
  0.0951361371292365,}
		} };

		static constexpr bool IsAdaptive = true;
		static constexpr size_t CorrectionMethodsCount = 1;
		static constexpr std::array<size_t, CorrectionMethodsCount> CorrectionMethodsAccuracyOrders = { 7 };
		static constexpr std::array<std::array<double, StepCount>, CorrectionMethodsCount> ButcherTableauErrorRow =
		{ {
			{
				-7.259091782802626e-5, 0, 0, 0, 0,
  -0.0010728916072503584,
  0.0002666668345794398,
  2.091533979096395,
  0.3213186752428666,
  -0.921013671395284,
  1.4012004409899175,
  0.0951361371292365,
  -2.9872967453726327,
			}
		} };

		static constexpr bool IsDenseOutputSupported = false;
	};
}