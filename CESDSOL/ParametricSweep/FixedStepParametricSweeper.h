#pragma once

#include "ParametricSweep/ParametricSweeper.h"
#include "Utils/Aliases.h"
#include "Utils/EventExecutor.h"
#include "Utils/Logger.h"
#include "Utils/Utils.h"

namespace CESDSOL
{
	enum class FixedStepParametricSweeperEvent
	{
		StartSweep,
		StartSolution,
		SuccessfulSolution,
		FailedSolution,
		FinishSweep,
		Count
	};

	template<typename ProblemType>
	class FixedStepParametricSweeper final
		: public ParametricSweeper<ProblemType>
		, public EventExecutor<FixedStepParametricSweeper<ProblemType>, ProblemType, FixedStepParametricSweeperEvent>
	{
	public:
		using ParametricSweeper<ProblemType>::ValueType;

		template<typename SolverType>
		FixedStepParametricSweeper
			( sptr<ProblemType> aProblem
			, uptr<SolverType> aSolver
			, size_t aParameterIndex = 0
			, ValueType aInitialValue = 0
			, ValueType aFinalValue = 0
			, ValueType aStep = 0
			) noexcept
			: ParametricSweeper<ProblemType>(std::move(aProblem), std::move(aSolver))
			, initialValue(aInitialValue)
			, finalValue(aFinalValue)
			, step(aStep)
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
			ValueType currentStep = finalValue < initialValue ? -step : step;

			const std::string& problemName = this->problem->GetDescriptor().GetProblemName();
			const std::string& parameterName = this->problem->GetDescriptor().GetParameterName(parameterIndex);
			
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
				Format("Starting parametric sweep over {} parameter of problem {}.", parameterName, problemName));
			
			this->ApplyActions(FixedStepParametricSweeperEvent::StartSweep, *this->problem);
			ValueType parameter = initialValue;
			ValueType oldStep = currentStep;
			Array<ValueType> previousSolution, tmp;
			if (interpolateInitialGuess)
			{
				previousSolution = this->problem->GetVariables().Flatten();
			}
			size_t stepIndex = 0;
			while (true)
			{
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
					Format("Starting solution with {} = {}.", parameterName, parameter));
				
				this->problem->SetParameter(parameterIndex, parameter);
				this->ApplyActions(FixedStepParametricSweeperEvent::StartSolution, *this->problem);
				const auto result = this->solver->Solve(*this->problem);
				if (!result->success)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::ParametricSweeper,
						Format("Parametric sweep stopped at {} = {} due to solver failure.", parameterName, parameter));
					
					this->ApplyActions(FixedStepParametricSweeperEvent::FailedSolution, *this->problem);
					break;
				}
				this->ApplyActions(FixedStepParametricSweeperEvent::SuccessfulSolution, *this->problem);
				stepIndex++;
				if (finalValue == parameter)
				{
					break;
				}
				const bool isFinal = (finalValue - parameter) * (finalValue - parameter - currentStep) < 0;
				oldStep = currentStep;
				if (isFinal)
				{
					currentStep = finalValue - parameter;
				}
				parameter += currentStep;
				if (interpolateInitialGuess)
				{
					tmp = previousSolution;
					previousSolution = this->problem->GetVariables().Flatten();
					AXPBY(-currentStep / oldStep, tmp, (oldStep + currentStep) / oldStep, this->problem->GetVariables().Flatten());
					this->problem->SetVariablesUpdated();
				}
			}
			this->ApplyActions(FixedStepParametricSweeperEvent::FinishSweep, *this->problem);
			return std::make_unique<OutputInfo>(finalValue == parameter, parameter);
		}
		
	private:
		size_t parameterIndex;

		MakeProperty(initialValue, InitialValue, ValueType, 0);
		MakeProperty(finalValue, FinalValue, ValueType, 0);
		MakeProperty(step, Step, ValueType, 0);
		MakeProperty(interpolateInitialGuess, InterpolateInitialGuess, bool, true);
	};
}