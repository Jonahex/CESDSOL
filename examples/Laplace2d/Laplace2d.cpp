#include "CESDSOL.h"

using namespace CESDSOL;

int main()
{
	// Autogen start
	auto descriptor = StationaryProblemDescriptor<2, double, double>(GridDescriptor<2, double>(5), Array<Array<std::array<size_t, 2>>>{ { {0, 2}, { 2, 0 }}}, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	descriptor.SetProblemName("Laplace2d");
	descriptor.SetParameterName(0, "k");
	descriptor.SetVariableName(0, "u");
	descriptor.SetContinuousEquation(0, 0, [](const auto& l, const auto& g) {return l.DerivativeValues[0][0] + l.DerivativeValues[0][1]; });
	descriptor.SetContinuousEquation(0, 1, [](const auto& l, const auto& g) {return l.FieldValues[0]; });
	descriptor.SetContinuousEquation(0, 2, [](const auto& l, const auto& g) {return l.FieldValues[0]; });
	descriptor.SetContinuousEquation(0, 3, [](const auto& l, const auto& g) {return -(g.Parameters[0] * (2 - l.Point[0]) * l.Point[0]) + l.FieldValues[0]; });
	descriptor.SetContinuousEquation(0, 4, [](const auto& l, const auto& g) {return l.FieldValues[0]; });
	descriptor.SetJacobianComponent(0, 0, 1, 0, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 2, 0, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 1, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 2, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 3, [](const auto& l, const auto& g) {return 1; });
	descriptor.SetJacobianComponent(0, 0, 0, 4, [](const auto& l, const auto& g) {return 1; });
	auto grid = std::make_shared<DirectProductGrid<2, double>>(SingleDimensionalGrid<double>({ MakeUniformRange<double>(0, 2, 100), std::nullopt }), SingleDimensionalGrid<double>({ MakeUniformRange<double>(0, 1, 100), std::nullopt }));
	std::unique_ptr<Discretization<2>> discretization = std::make_unique<StructuredFiniteDifferenceDiscretization<2>>(5);
	auto problem = descriptor.MakeProblem(grid, std::move(discretization));
	// Autogen end

	problem->SetVariables(Vector<double>(0., problem->DOFCount()));

	auto pardiso = std::make_unique<MKL::PARDISO<double>>();
	auto gss = MakeLineSearcher<GoldenSectionSearch>(*problem);
	auto newton = MakeNonlinearSolver<ModifiedNewton>(*problem, std::move(pardiso), std::move(gss));

	auto sweeper = FixedStepParametricSweeper(problem, std::move(newton), 0, 1., 5, 1.);
	sweeper.AddAction(FixedStepParametricSweeperEvent::SuccessfulSolution, [](auto& problem) {
		Serializer::Save(problem, std::format("{}Laplace2d/output/{}.dat", ExamplesPath, problem.MakeStateName())); });
	sweeper.Sweep();
}