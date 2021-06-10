#include "CESDSOL.h"

using namespace CESDSOL;

int main()
{
	// Autogen start
	auto descriptor = ExplicitTransientProblemDescriptor<1, double, double>(GridDescriptor<1, double>(3), Array<Array<std::array<size_t, 1>>>{ { {2}}, {}, { {1} }, { {1} }, { {1} }, { {1} }}, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	descriptor.SetProblemName("FermionSineGordonKinkCollision");
	descriptor.SetParameterName(0, "g");
	descriptor.SetVariableName(0, "f");
	descriptor.SetVariableName(1, "f$t1");
	descriptor.SetVariableName(2, "u");
	descriptor.SetVariableName(3, "v");
	descriptor.SetVariableName(4, "p");
	descriptor.SetVariableName(5, "q");
	descriptor.SetContinuousEquation(0, 0, [](const auto& l, const auto& g) {return l.FieldValues[1]; });
	descriptor.SetContinuousEquation(1, 0, [](const auto& l, const auto& g) {return -sin(l.FieldValues[0]) + 2 * g.Parameters[0] * (l.FieldValues[4] * l.FieldValues[5] + l.FieldValues[2] * l.FieldValues[3]) + l.DerivativeValues[0][0]; });
	descriptor.SetContinuousEquation(2, 0, [](const auto& l, const auto& g) {return g.Parameters[0] * l.FieldValues[0] * l.FieldValues[5] - l.DerivativeValues[5][0]; });
	descriptor.SetContinuousEquation(3, 0, [](const auto& l, const auto& g) {return g.Parameters[0] * l.FieldValues[0] * l.FieldValues[4] + l.DerivativeValues[4][0]; });
	descriptor.SetContinuousEquation(4, 0, [](const auto& l, const auto& g) {return -(g.Parameters[0] * l.FieldValues[0] * l.FieldValues[3]) + l.DerivativeValues[3][0]; });
	descriptor.SetContinuousEquation(5, 0, [](const auto& l, const auto& g) {return -(g.Parameters[0] * l.FieldValues[0] * l.FieldValues[2]) - l.DerivativeValues[2][0]; });
	auto grid = std::make_shared<DirectProductGrid<1, double>>(SingleDimensionalGrid<double>({ MakeUniformRange<double>(-100, 100, 2000), std::nullopt }));
	std::unique_ptr<Discretization<1>> discretization = std::make_unique<StructuredFiniteDifferenceDiscretization<1>>(5);
	auto problem = descriptor.MakeProblem(grid, std::move(discretization));
	problem->AddLocalOutputExpression([](const auto& l, const auto& g) { return sqrt(pow(l.FieldValues[4], 2) + pow(l.FieldValues[5], 2) + pow(l.FieldValues[2], 2) + pow(l.FieldValues[3], 2)); }, "norm");
	// Autogen end

	// Initial value.
	const size_t gridSize = 2001;
	const double omega = 1;
	const double sigma = 1;
	const double x0 = -20;
	const double time = 100;

	for (size_t gridIndex = 0; gridIndex < gridSize; ++gridIndex)
	{
		const auto A = 1 / (sqrt(2 * sqrt(Pi * sigma)));
		const auto x = grid->GetCoordinates(gridIndex)[0];
		const auto delta = x - x0;
		const auto a = A * exp(-(delta * delta) / (2 * sigma)) * cos(omega * delta);
		const auto b = A * exp(-(delta * delta) / (2 * sigma)) * sin(omega * delta);

		problem->SetVariable(0, gridIndex, 4 * atan(exp(x)));
		problem->SetVariable(1, gridIndex, 0);
		problem->SetVariable(2, gridIndex, a);
		problem->SetVariable(3, gridIndex, -b);
		problem->SetVariable(4, gridIndex, b);
		problem->SetVariable(5, gridIndex, a);
	}
	problem->SetVariablesUpdated();
	problem->SetParameter(0, 1);

	auto rk = RungeKuttaSolver<Verner87, double, double, true, true>();
	rk.SetAbsoluteTolerance(1e-8);
	rk.SetRelativeTolerance(1e-8);
	rk.SetDenseOutputStep(1.);
	rk.Solve(0, time, *problem);

	Serializer::Save(*problem, std::format("{}FermionSineGordonKinkCollision/output/{}.dat", ExamplesPath, problem->MakeStateName()));
}