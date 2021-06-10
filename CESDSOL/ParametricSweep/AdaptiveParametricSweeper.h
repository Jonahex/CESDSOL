#pragma once

#include "ParametricSweep/ParametricSweeper.h"
#include "Utils/EventExecutor.h"

namespace CESDSOL
{
	enum class AdaptiveParametricSweeperEvent
	{
		StartSweep,
		StartSolution,
		SuccessfulSolution,
		FailedSolution,
		FinishSweep,
		StartBranchChange,
		FailedBranchChangeAttempt,
		Count
	};

	template<typename ProblemType>
	class AdaptiveParametricSweeper final
		: public ParametricSweeper<ProblemType>
		, public EventExecutor<AdaptiveParametricSweeper<ProblemType>, ProblemType, AdaptiveParametricSweeperEvent>
	{
	public:
		using ParametricSweeper<ProblemType>::ValueType;

		template<typename SolverType>
		AdaptiveParametricSweeper
			( sptr<ProblemType> aProblem
			, uptr<SolverType> aSolver
			, size_t aParameterIndex = 0
			, ValueType aInitialValue = 0
			, ValueType aFinalValue = 0
			, ValueType aInitialStep = 0
			) noexcept
			: ParametricSweeper<ProblemType>(std::move(aProblem), std::move(aSolver))
			, initialValue(aInitialValue)
			, finalValue(aFinalValue)
			, initialStep(aInitialStep)
		{
			SetParameterIndex(aParameterIndex);
		}

		void SetParameterIndex(size_t aParameterIndex) noexcept
		{
			AssertE(aParameterIndex < this->problem->ParameterCount(), MessageTag::ParametricSweeper,
				Format("Parameter index %zu exceeds problem parameter count %zu.", aParameterIndex, this->problem->ParameterCount()));
			parameterIndex = aParameterIndex;
		}

		[[nodiscard]] size_t GetParameterIndex() const noexcept
		{
			return parameterIndex;
		}

		struct OutputInfo : ParametricSweeper<ProblemType>::OutputInfo
		{
			OutputInfo(bool aSuccess, ValueType aFinalValue)
				: ParametricSweeper<ProblemType>::OutputInfo(aSuccess)
				, finalValue(aFinalValue)
			{}

			ValueType finalValue = 0;
		};

		uptr<typename ParametricSweeper<ProblemType>::OutputInfo> Sweep() const noexcept override
		{
			const std::string& problemName = this->problem->GetDescriptor().GetProblemName();
			const std::string& parameterName = this->problem->GetDescriptor().GetParameterName(parameterIndex);

			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
				Format("Starting parametric sweep over {} parameter of problem {}.", parameterName, problemName));

			this->ApplyActions(AdaptiveParametricSweeperEvent::StartSweep, *this->problem);

			ValueType currentStep = initialStep, oldStep = currentStep;
			ValueType parameter = initialValue, previousParameter = initialValue;
			Vector<ValueType> previousSolution, tmp;
			size_t solutionIndex = 0;
			size_t branch = 0;
			bool isChangingBranch = false;
			ValueType changeBranchStep = 1;
			size_t changeBranchTrial = 0;

			if (interpolateInitialGuess)
			{
				previousSolution = this->problem->GetVariables().Flatten();
				tmp = previousSolution;
			}

			while (true)
			{
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
					Format("Starting solution with {} = {} on branch {}.", parameterName, parameter, branch + 1));

				this->problem->SetParameter(parameterIndex, parameter);
				this->ApplyActions(AdaptiveParametricSweeperEvent::StartSolution, *this->problem);
				const auto result = this->solver->Solve(*this->problem);
				if (!result->success)
				{
					if (isChangingBranch && changeBranchTrial > 0 && changeBranchTrial <= maxChangeBranchTrials)
					{
						Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::ParametricSweeper,
							Format("Parametric sweep failed to change branch at {} = {}. Trying again...",
								parameterName, parameter));
						this->ApplyActions(AdaptiveParametricSweeperEvent::FailedBranchChangeAttempt, *this->problem);

						changeBranchStep *= growthFactor;
						changeBranchTrial++;
						this->problem->SetVariables(previousSolution);
						previousSolution = tmp;
						if (interpolateInitialGuess)
						{
							AXPBY(-changeBranchStep, previousSolution, 1 + changeBranchStep, this->problem->GetVariables().Flatten());
							this->problem->SetVariablesUpdated();
						}
						oldStep = currentStep;
						currentStep *= growthFactor;
						parameter = previousParameter + currentStep;
					}
					else
					{
						this->ApplyActions(AdaptiveParametricSweeperEvent::FailedSolution, *this->problem);

						currentStep /= shrinkFactor;

						if (std::abs(currentStep) < minStep)
						{
							if (isChangingBranch || !tryChangeBranch)
							{
								Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
									Format("Parametric sweep stopped at {} = {} on branch {} due to step underflow.", 
										parameterName, parameter, branch + 1));
								break;
							}
							else
							{
								Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::ParametricSweeper,
									Format("Parametric sweep is trying to change branch at {} = {} due to step underflow.",
										parameterName, parameter));
								this->ApplyActions(AdaptiveParametricSweeperEvent::StartBranchChange, *this->problem);

								changeBranchTrial++;

								isChangingBranch = true;
								this->problem->SetVariables(previousSolution);
								previousSolution = tmp;
								if (interpolateInitialGuess)
								{
									AXPBY(-changeBranchStep, previousSolution, 1 + changeBranchStep, this->problem->GetVariables().Flatten());
									this->problem->SetVariablesUpdated();
								}
								oldStep = currentStep;
								currentStep *= -1;
								parameter = previousParameter + currentStep;
							}
						}
						else
						{
							Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::ParametricSweeper,
								Format("Parametric sweep is decreasing step at {} = {} on branch {} due to solver failure.",
									parameterName, parameter, branch + 1));

							this->problem->SetVariables(previousSolution);
							parameter = previousParameter;
						}
					}
				}
				else
				{
					solutionIndex++;
					if (isChangingBranch)
					{
						branch++;
					}
					oldStep = currentStep;
					currentStep = currentStep / std::abs(currentStep) * std::min(std::abs(currentStep * growthFactor), maxStep);
					this->ApplyActions(AdaptiveParametricSweeperEvent::SuccessfulSolution, *this->problem);

					if (limitBranchCount && branch > maxBranchCount)
					{
						Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
							Format("Parametric sweep finished not reaching target paramter value at {} = {} on branch {} due to reaching maximum branch count.",
								parameterName, parameter, branch + 1));
						break;
					}
					if (limitSolutionCount && solutionIndex > maxSolutionCount)
					{
						Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
							Format("Parametric sweep finished not reaching target parameter value at {} = {} on branch {} due to reaching maximum solution count.",
								parameterName, parameter, branch + 1));
						break;
					}

					tmp = previousSolution;
					previousSolution = this->problem->GetVariables().Flatten();

					if (isChangingBranch)
					{
						changeBranchStep = 1;
						changeBranchTrial = 0;
						isChangingBranch = false;
					}
					else if (interpolateInitialGuess)
					{
						AXPBY(-currentStep / oldStep, tmp, (oldStep + currentStep) / oldStep, this->problem->GetVariables().Flatten());
						this->problem->SetVariablesUpdated();
					}
				}

				if (parameter == finalValue)
				{
					break;
				}

				if (!isChangingBranch)
				{
					previousParameter = parameter;
					const bool isFinal = (finalValue - parameter) * (finalValue - parameter - currentStep) < 0;
					if (isFinal)
					{
						currentStep = finalValue - parameter;
					}
					parameter += currentStep;
				}
			}
			this->ApplyActions(AdaptiveParametricSweeperEvent::FinishSweep, *this->problem);
			return std::make_unique<OutputInfo>(finalValue == parameter, parameter);
		}

	private:
		size_t parameterIndex;

		MakeProperty(initialValue, InitialValue, ValueType, 0);
		MakeProperty(finalValue, FinalValue, ValueType, 0);
		MakeProperty(initialStep, InitialStep, ValueType, 0.01);
		MakeProperty(minStep, MinStep, ValueType, 1e-6);
		MakeProperty(maxStep, MaxStep, ValueType, 1e-1);
		MakeProperty(interpolateInitialGuess, InterpolateInitialGuess, bool, true);
		MakeProperty(growthFactor, GrowthFactor, ValueType, 1.1);
		MakeProperty(shrinkFactor, ShrinkFactor, ValueType, 1.5);
		MakeProperty(tryChangeBranch, TryChangeBranch, bool, true);
		MakeProperty(maxChangeBranchTrials, MaxChangeBranchTrials, size_t, 5);
		MakeProperty(limitBranchCount, LimitBranchCount, bool, true);
		MakeProperty(maxBranchCount, MaxBranchCount, size_t, 5);
		MakeProperty(limitSolutionCount, LimitSolutionCount, bool, true);
		MakeProperty(maxSolutionCount, MaxSolutionCount, size_t, 1000);
	};
}