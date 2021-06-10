#include "CESDSOL.h"

using namespace CESDSOL;

int main()
{
	// Autogen start
	auto descriptor = StationaryProblemDescriptor<1, double, double>(GridDescriptor<1, double>(3), Array<Array<std::array<size_t, 1>>>{ { {2}}}, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	descriptor.SetProblemName("StaticSineGordon");
	descriptor.SetParameterName(0, "g");
	descriptor.SetVariableName(0, "phi");
	descriptor.SetContinuousEquation(0, 0, [](const auto& l, const auto& g) {return -(g.Parameters[0] * sin(l.FieldValues[0])) + l.DerivativeValues[0][0]; });
	descriptor.SetContinuousEquation(0, 1, [](const auto& l, const auto& g) {return -Pi + l.FieldValues[0]; });
	descriptor.SetContinuousEquation(0, 2, [](const auto& l, const auto& g) {return -2 * Pi + l.FieldValues[0]; });
	descriptor.SetJacobianComponent(0, 0, 0, 0, [](const auto& l, const auto& g) {return -(g.Parameters[0] * cos(l.FieldValues[0])); });
	descriptor.SetJacobianComponent(0, 0, 1, 0, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 1, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 2, [](const auto& l, const auto& g) {return 1; });
	auto grid = std::make_shared<DirectProductGrid<1, double>>(SingleDimensionalGrid<double>({ MakeUniformRange<double>(0, 10, 1000), std::nullopt }));
	std::unique_ptr<Discretization<1>> discretization = std::make_unique<StructuredFiniteDifferenceDiscretization<1>>(5);
	auto problem = descriptor.MakeProblem(grid, std::move(discretization));
	// Autogen end

	problem->SetVariables(Vector<double>(2 * Pi, problem->DOFCount()));

	auto pardiso = std::make_unique<MKL::PARDISO<double>>();
	auto gss = MakeLineSearcher<GoldenSectionSearch>(*problem);
	auto newton = MakeNonlinearSolver<ModifiedNewton>(*problem, std::move(pardiso), std::move(gss));

	auto sweeper = FixedStepParametricSweeper(problem, std::move(newton), 0, 1., 10, 1.);
	sweeper.AddAction(FixedStepParametricSweeperEvent::SuccessfulSolution, [](auto& problem) {
		Serializer::Save(problem, std::format("{}StaticSineGordon/output/{}.dat", ExamplesPath, problem.MakeStateName())); });
	sweeper.Sweep();
}