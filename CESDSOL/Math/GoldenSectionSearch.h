#pragma once

#include "Math/LineSearcher.h"
#include "Math/VectorOperations.h"
#include "Utils/Utils.h"

#include <chrono>

namespace CESDSOL
{
	enum class GSSExitConditions : uint32_t
	{
		MeritGoalReached = 1 << 0,
		IterationCount = 1 << 1,
		SolutionStagnation = 1 << 2,
		MeritStagnation = 1 << 3,
		Everything = (1 << 4) - 1
	};
	MakeFlag(GSSExitConditions)
	
	template<typename ProblemType>
	class GoldenSectionSearch final
		: public LineSearcher<ProblemType>
	{
	private:
		using LineSearcher<ProblemType>::ValueType;
		
		static ValueType GetMerit(ProblemType& problem, const Vector<ValueType>& previousSolution,
			const Vector<ValueType>& shift, double currentMultiplier) noexcept
		{
			problem.SetVariables(previousSolution);
			AXPY(currentMultiplier, shift, problem.GetVariables().Flatten());
			return problem.GetMerit();
		}

	public:
		MakeProperty(exitConditions, ExitConditions, GSSExitConditions, GSSExitConditions::Everything)
		MakeProperty(left, LeftBoundary, double, 0.)
		MakeProperty(right, RightBoundary, double, 1.)
		MakeProperty(solutionTolerance, SolutionTolerance, double, 1e-8)
		MakeProperty(meritTolerance, MeritTolerance, double, 1e-8)
		MakeProperty(meritGoal, MeritGoal, double, 1e-8)
		MakeProperty(iterationLimit, IterationLimit, size_t, 100)

		struct OutputInfo final : LineSearcher<ProblemType>::OutputInfo
		{
			OutputInfo(bool aSuccess, double aFinalMerit, size_t aIterationCount)
				: LineSearcher<ProblemType>::OutputInfo(aSuccess)
				, finalMerit(aFinalMerit)
				, iterationCount(aIterationCount)
			{}
			
			double finalMerit;
			size_t iterationCount;
		};

		GoldenSectionSearch() = default;

		GoldenSectionSearch(double aLeft, double aRight) noexcept
			: left(aLeft)
			, right(aRight)
		{}

		uptr<typename LineSearcher<ProblemType>::OutputInfo>
			Solve(ProblemType& problem, const Vector<ValueType>& shift) const noexcept override
		{
			Vector<ValueType> previousSolution = problem.GetVariables().Flatten();
			
			std::chrono::high_resolution_clock clock;
			const auto solutionStartTime = clock.now();

			Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::LineSearcher,
				"Starting line search using golden section method.");
			auto a = left, b = right, c = b - (b - a) / GoldenRatio, d = a + (b - a) / GoldenRatio;
			auto fc = GetMerit(problem, previousSolution, shift, c),
				fd = GetMerit(problem, previousSolution, shift, d);
			double fCurrent;
			size_t iterationNumber = 0;
			while (true)
			{
				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::LineSearcher,
					Format("Starting golden section search iteration {}:", iterationNumber));
				const auto iterationStartTime = clock.now();
				if (fc < fd)
				{
					b = d;
					d = c;
					fd = fc;
					c = b - (b - a) / GoldenRatio;
					fc = GetMerit(problem, previousSolution, shift, c);
					fCurrent = fc;
				}
				else
				{
					a = c;
					c = d;
					fc = fd;
					d = a + (b - a) / GoldenRatio;
					fd = GetMerit(problem, previousSolution, shift, d);
					fCurrent = fd;
				};

				if (exitConditions & GSSExitConditions::MeritGoalReached && fCurrent < meritGoal)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LineSearcher,
						Format("Golden section search converged after {} iterations in {}.", iterationNumber + 1, clock.now() - solutionStartTime));
					break;
				};
				if (exitConditions & GSSExitConditions::SolutionStagnation && b - a < solutionTolerance)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LineSearcher,
						Format("Stopping golden section search after {} iterations due to search range {} decreased below tolerance.", iterationNumber + 1, b - a));
					break;
				};
				if (exitConditions & GSSExitConditions::MeritStagnation && std::abs(fd - fc) < meritTolerance)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LineSearcher,
						Format("Stopping golden section search after {} iterations due to merit value change {} decreased below tolerance.", iterationNumber + 1, std::abs(fd - fc)));
					break;
				};

				++iterationNumber;
				if (exitConditions & GSSExitConditions::IterationCount && iterationNumber > iterationLimit)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LineSearcher,
						Format("Golden section search failed to find satisfying solution: maximum number of iterations {} is exceeded.", iterationLimit));
					break;
				};

				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::LineSearcher,
					Format("Finishing golden section search iteration with merit value {} in {}. Shift factor is {}.\n", fCurrent, clock.now() - iterationStartTime, (c + d) / 2));
			};

			return std::make_unique<OutputInfo>(true, fCurrent, iterationNumber);
		}
	};
}