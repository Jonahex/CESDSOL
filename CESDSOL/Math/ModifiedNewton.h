#pragma once

#include "Math/LinearSolver.h"
#include "Math/LineSearcher.h"
#include "Math/NonlinearSolver.h"
#include "Math/VectorOperations.h"
#include "Utils/Aliases.h"
#include "Utils/Logger.h"
#include "Utils/Utils.h"

#include <chrono>
#include <optional>

namespace CESDSOL
{
	enum class MNExitConditions : uint32_t
	{
		MeritGoalReached = 1 << 0,
		IterationCount = 1 << 1,
		MeritOverflow = 1 << 2,
		SolutionStagnation = 1 << 3,
		MeritStagnation = 1 << 4,
		MeritIncrease = 1 << 5,
		Everything = (1 << 6) - 1
	};
	MakeFlag(MNExitConditions)
	
	template<typename ProblemType>
	class ModifiedNewton final
		: public NonlinearSolver<ProblemType>
	{
	public:
		using CurrentLinearSolver = LinearSolver<typename ProblemType::JacobianMatrixType, typename ProblemType::VectorType>;
		using CurrentLineSearcher = LineSearcher<ProblemType>;
	
	private:
		uptr<CurrentLinearSolver> linearSolver;
		uptr<CurrentLineSearcher> lineSearcher;
		
	public:		
		MakeProperty(exitConditions, ExitConditions, MNExitConditions, MNExitConditions::MeritGoalReached
			| MNExitConditions::SolutionStagnation
			| MNExitConditions::MeritStagnation
			| MNExitConditions::IterationCount
			| MNExitConditions::MeritIncrease)
		MakeProperty(meritGoal, MeritGoal, double, 1e-8)
		MakeProperty(iterationLimit, IterationLimit, size_t, 100)
		MakeProperty(maxMerit, MaximumMerit, double, 1e10)
		MakeProperty(solutionTolerance, SolutionTolerance, double, 1e-10)
		MakeProperty(meritTolerance, MeritTolerance, double, 1e-10)
		MakeProperty(meritIncreaseFactor, MeritIncreaseFactor, double, 1)
		MakeProperty(dampring, Damping, double, 1)

		struct OutputInfo final : NonlinearSolver<ProblemType>::OutputInfo
		{
			OutputInfo(bool aSuccess, double aFinalMerit, size_t aIterationCount)
				: NonlinearSolver<ProblemType>::OutputInfo(aSuccess)
				, finalMerit(aFinalMerit)
				, iterationCount(aIterationCount)
			{}

			double finalMerit;
			size_t iterationCount;
		};

		ModifiedNewton(uptr<CurrentLinearSolver> aLinearSolver, uptr<CurrentLineSearcher> aLineSearcher)
			: linearSolver(std::move(aLinearSolver))
			, lineSearcher(std::move(aLineSearcher))
		{}

		uptr<typename NonlinearSolver<ProblemType>::OutputInfo>
			Solve(ProblemType& problem) const noexcept override
		{
			using ValueType = typename ProblemType::FieldType;

			auto tmp = Vector<ValueType>(problem.DOFCount());
			ValueType oldMerit = 0, oldSolutionNorm;
			Vector<ValueType> oldSolution;
			std::chrono::high_resolution_clock clock;
			size_t iterationCount = 0;

			Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::NonlinearSolver, 
				"Starting solution of nonlinear equation system using modified Newton method.\n");
			const auto solutionStartTime = clock.now();
			while (true)
			{
				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::NonlinearSolver,
					Format("Starting modified Newton iteration {}:", iterationCount + 1));
				const auto iterationStartTime = clock.now();
				if (!linearSolver->Solve(problem.GetJacobian(), problem.GetEquations().Flatten(), tmp))
				{
					Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
						Format("Stopping modified Newton solution due to linear solver failure after {} iterations.", iterationCount + 1));
					return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
				}

				Scale(-1, tmp);
				const auto result = lineSearcher->Solve(problem, tmp);
				if (!result->success)
				{
					Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
						Format("Stopping modified Newton solution due to line searcher failure after {} iterations.", iterationCount + 1));
					return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
				}

				const auto merit = problem.GetMerit();
				const auto solutionNorm = problem.CalculateSolutionNorm();
				if (exitConditions & MNExitConditions::MeritGoalReached && merit < meritGoal)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::NonlinearSolver,
						Format("Modified Newton solver successfully converged after {} iterations in {}.\n", iterationCount + 1, clock.now() - solutionStartTime));
					return std::make_unique<OutputInfo>(true, oldMerit, iterationCount);
				}
				if (exitConditions & MNExitConditions::MeritOverflow && merit > maxMerit)
				{
					Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
						Format("Stopping modified Newton solution after {} iterations due to merit overflow with merit value {}.\n", iterationCount + 1, merit));
					return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
				}
				if (iterationCount > 0)
				{
					if (exitConditions & MNExitConditions::MeritIncrease && merit > meritIncreaseFactor * oldMerit)
					{
						Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
							Format("Stopping modified Newton solution on {} iteration due to merit increased {} times since last iteration.\n", iterationCount + 1, merit / oldMerit));
						return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
					}
					if (exitConditions & MNExitConditions::MeritStagnation && std::abs(merit - oldMerit) < meritTolerance)
					{
						Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
							Format("Stopping modified Newton solution after {} iterations due to merit value change {} decreased below tolerance.\n", iterationCount + 1, std::abs(merit - oldMerit)));
						return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
					}
					if (exitConditions & MNExitConditions::SolutionStagnation && std::abs(solutionNorm - oldSolutionNorm) < solutionTolerance)
					{
						Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
							Format("Stopping modified Newton solution after {} iterations due to solution norm change {} decreased below tolerance.\n", iterationCount + 1, std::abs(solutionNorm - oldSolutionNorm)));
						return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
					}
				}
				oldMerit = merit;
				oldSolutionNorm = solutionNorm;

				++iterationCount;
				if (exitConditions & MNExitConditions::IterationCount && iterationCount > iterationLimit)
				{
					Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::NonlinearSolver,
						Format("Stopping modified Newton solution on iteration {} due to hitting iteration limit with merit value {}.\n", iterationCount + 1, merit));
					return std::make_unique<OutputInfo>(false, oldMerit, iterationCount);
				}

				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::NonlinearSolver,
					Format("Finishing modified Newton iteration with merit value {} in {}.\n", merit, clock.now() - iterationStartTime));
			}
		}
	};
}