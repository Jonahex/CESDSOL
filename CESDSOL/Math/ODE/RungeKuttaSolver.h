#pragma once

#include "Math/LinearAlgebra.h"

namespace CESDSOL
{
	template<typename ParameterType = double, typename ValueType = double, bool UseAdaptive = true, bool UseDenseOutput = true>
	class RKOptions {};

	template<typename ParameterType, typename ValueType>
	class RKOptions<ParameterType, ValueType, true, false>
	{
		MakeProperty(absoluteTolerance, AbsoluteTolerance, ParameterType, 1e-10);
		MakeProperty(relativeTolerance, RelativeTolerance, ParameterType, 1e-10);
		MakeProperty(minStep, MinStep, ParameterType, 1e-12);
		MakeProperty(maxStepScale, MaxStepScale, ParameterType, 6.0);
		MakeProperty(minStepScale, MinStepScale, ParameterType, 0.333);
		MakeProperty(maxError, MaxError, ParameterType, 1.);
		MakeProperty(stepScaleFactor, StepScaleFactor, ParameterType, 0.9);
	};

	template<typename ParameterType, typename ValueType, bool UseAdaptive>
	class RKOptions<ParameterType, ValueType, UseAdaptive, true> : public RKOptions<ParameterType, ValueType, UseAdaptive, false>
	{
		MakeProperty(denseOutputStep, DenseOutputStep, ParameterType, 1e-3);
	};

	enum class RKExitConditions : uint32_t
	{
		StepUnderflow = 1 << 0,
		StepCountLimitReached = 1 << 1,
		SolutionNormOverflow = 1 << 2,
		Everything = (1 << 3) - 1
	};
	MakeFlag(RKExitConditions)
	
	template<typename MethodDescriptor, typename ParameterType = double, typename ValueType = double, bool UseAdaptive = true, bool UseDenseOutput = true>
	class RungeKuttaSolver : public RKOptions<ParameterType, ValueType, UseAdaptive && MethodDescriptor::IsAdaptive, UseDenseOutput && MethodDescriptor::IsDenseOutputSupported>
	{
	private:
		static constexpr bool ActualUseAdaptive = UseAdaptive && MethodDescriptor::IsAdaptive;
		static constexpr bool ActualUseDenseOutput = UseDenseOutput && MethodDescriptor::IsDenseOutputSupported;

	public:
		MakeProperty(exitConditions, ExitConditions, RKExitConditions, RKExitConditions::Everything)
		MakeProperty(initialStep, InitialStep, ParameterType, 1e-3);
		MakeProperty(stepCountLimit, StepCountLimit, size_t, 50000);
		MakeProperty(maxSolutionNorm, MaxSolutionNorm, ValueType, 1e20);

		struct OutputInfo
		{
			enum class OutputReason
			{
				Success,
				StepUnderflow,
				StepCountLimitReached,
				SolutionNormOverflow
			};

			bool success;
			OutputReason reason;
			size_t stepCount;
		};

	private:
		struct SolutionTemporaries
		{
			size_t stepCount = 0;

			ParameterType currentX;
			ParameterType previousX;
			ParameterType currentStep;
			ParameterType previousStep;
			Vector<ValueType> currentY;
			Vector<ValueType> previousY;
			ValueType currentError = 1;
			ValueType previousError;
			Vector<ValueType> currentEquations;
			Vector<ValueType> previousEquations;

			std::array<Vector<ValueType>, MethodDescriptor::StepCount> stepArrays;

			SolutionTemporaries(ParameterType firstX, const Vector<ValueType>& firstY, ParameterType aInitialStep)
				: currentX(firstX)
				, currentStep(aInitialStep)
				, currentY(firstY)
				, previousEquations(firstY.size())
			{
				for (auto& stepArray : stepArrays)
				{
					stepArray = Vector<ValueType>(currentY.size());
				}
			}
		};
		
		template<typename ProblemType, bool UseDenseOutput>
		class DenseOutputImpl
		{
		public:
			DenseOutputImpl(const SolutionTemporaries& tmp) {}
		};

		template<typename ProblemType>
		class DenseOutputImpl<ProblemType, true>
		{
		private:
			[[nodiscard]] static auto ComputePolynomialCoefficients(ParameterType theta) noexcept
			{
				std::array<ParameterType, MethodDescriptor::StepCount + MethodDescriptor::DenseOutputStepCount + 1> result{};
				ParameterType mult = 1;
				for (size_t power = 0; power <= MethodDescriptor::InterpolationOrder; ++power)
				{
					for (size_t stepIndex = 0; stepIndex < result.size(); ++stepIndex)
					{
						result[stepIndex] += mult * MethodDescriptor::DenseOutputCoefficients[stepIndex][power];
					}
					mult *= theta;
				}
				return result;
			}

			std::array<Vector<ValueType>, MethodDescriptor::DenseOutputStepCount> additionalSteps;
			size_t currentGridIndex = 0;

		public:
			DenseOutputImpl(const SolutionTemporaries& tmp)
			{
				for (auto& step : additionalSteps)
				{
					step = Vector<typename ProblemType::FieldType>(tmp.currentY.size());
				}
			}

			void Process(SolutionTemporaries& tmp, ParameterType denseOutputStep, ParameterType firstX,
				ProblemType& problem)
			{
				const size_t nextGridIndex = std::floor((tmp.currentX - firstX) / denseOutputStep);

				if (nextGridIndex > currentGridIndex)
				{
					if constexpr (MethodDescriptor::DenseOutputStepCount > 0)
					{
						for (size_t i = 0; i < MethodDescriptor::DenseOutputStepCount; i++)
						{
							for (size_t j = 0; j < problem.DOFCount(); j++)
							{
								ValueType sum = MethodDescriptor::DenseOutputButcherTableauMainPart[i][MethodDescriptor::StepCount] * tmp.currentEquations[j];
								for (size_t k = 0; k < MethodDescriptor::StepCount; k++)
								{
									sum += MethodDescriptor::DenseOutputButcherTableauMainPart[i][k] * tmp.stepArrays[k][j];
								}
								for (size_t k = 0; k < i; ++k)
								{
									sum += MethodDescriptor::DenseOutputButcherTableauMainPart[i][k + MethodDescriptor::StepCount + 1] * additionalSteps[k][j];
								}
								problem.GetVariables().Flatten()[j] = tmp.previousY[j] + tmp.previousStep * sum;
							}
							problem.SetTime(tmp.previousX + MethodDescriptor::DenseOutputButcherTableauFirstColumn[i] * tmp.previousStep);
							problem.SetVariablesUpdated();
							additionalSteps[i] = problem.GetEquations().Flatten();
						}
					}
					
					while (currentGridIndex < nextGridIndex)
					{
						currentGridIndex++;
						const auto theta = (firstX + denseOutputStep * currentGridIndex - tmp.previousX) / tmp.previousStep;
						const auto coefficients = ComputePolynomialCoefficients(theta);
						for (size_t equationIndex = 0; equationIndex < problem.DOFCount(); ++equationIndex)
						{
							ValueType sum = tmp.currentEquations[equationIndex] * coefficients[MethodDescriptor::StepCount];
							for (size_t stepIndex = 0; stepIndex < MethodDescriptor::StepCount; ++stepIndex)
							{
								sum += tmp.stepArrays[stepIndex][equationIndex] * coefficients[stepIndex];
							}
							for (size_t stepIndex = 0; stepIndex < MethodDescriptor::DenseOutputStepCount; ++stepIndex)
							{
								sum += additionalSteps[stepIndex][equationIndex] * coefficients[stepIndex + MethodDescriptor::StepCount + 1];
							}
							problem.GetVariables().Flatten()[equationIndex] = tmp.previousY[equationIndex] + tmp.previousStep * sum;
						}
						problem.SetTime(firstX + denseOutputStep * currentGridIndex);
						problem.CacheCurrent();
					}

					problem.SetTime(tmp.currentX);
				}					
			}
		};

		template<typename ProblemType>
		using DenseOutput = DenseOutputImpl<ProblemType, ActualUseDenseOutput>;

	public:
		template<typename ProblemType>
		[[nodiscard]] OutputInfo Solve(ParameterType firstX, ParameterType lastX, ProblemType& problem) const noexcept
		{			
			problem.SetTime(firstX);
			problem.CacheCurrent();

			SolutionTemporaries tmp{ firstX, problem.GetVariables().Flatten(), initialStep };
			DenseOutput<ProblemType> denseOutput{ tmp };

			std::chrono::high_resolution_clock clock;
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::TransientSolver,
				"Starting solution of transient problem using Runge-Kutta method.\n");
			const auto solutionStartTime = clock.now();

			tmp.currentEquations = problem.GetEquations().Flatten();
			while (tmp.currentX < lastX)
			{
				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::TransientSolver,
					Format("Starting calculation of new step at time {}:", tmp.currentX));
				const auto iterationStartTime = clock.now();
				tmp.stepArrays[0] = tmp.currentEquations;
				tmp.previousY = std::move(tmp.currentY);
				tmp.previousError = tmp.currentError;
				tmp.previousStep = tmp.currentStep;
				for (size_t i = 1; i <= MethodDescriptor::StepCount; i++)
				{
					for (size_t j = 0; j < problem.DOFCount(); j++)
					{
						ValueType sum = 0;
						for (size_t k = 0; k < i; k++)
						{
							sum += MethodDescriptor::ButcherTableauMainPart[i][k] * tmp.stepArrays[k][j];
						}
						problem.GetVariables().Flatten()[j] = tmp.previousY[j] + tmp.currentStep * sum;
					}
					problem.SetVariablesUpdated();
					if (i != MethodDescriptor::StepCount)
					{
						problem.SetTime(tmp.currentX + MethodDescriptor::ButcherTableauFirstColumn[i] * tmp.currentStep);
						tmp.stepArrays[i] = problem.GetEquations().Flatten();
					}
					else
					{
						problem.SetTime(tmp.currentX + tmp.currentStep);
						tmp.currentEquations = problem.GetEquations().Flatten();
					}
				}
				tmp.currentY = problem.GetVariables().Flatten();

				bool canContinue = true;
				if constexpr (ActualUseAdaptive)
				{
					std::array<ValueType, MethodDescriptor::CorrectionMethodsCount> error{};
					for (size_t j = 0; j < problem.DOFCount(); j++)
					{
						ValueType errorTemp = this->GetAbsoluteTolerance() + this->GetRelativeTolerance() * std::max(std::abs(tmp.currentY[j]), std::abs(tmp.previousY[j]));
						for (size_t i = 0; i < MethodDescriptor::CorrectionMethodsCount; i++)
						{
							ValueType sum = 0;
							for (size_t k = 0; k < MethodDescriptor::StepCount; k++)
							{
								sum += MethodDescriptor::ButcherTableauErrorRow[i][k] * tmp.stepArrays[k][j];
							}
							if constexpr (MethodDescriptor::ButcherTableauErrorRow.front().size() > MethodDescriptor::StepCount)
							{
								sum += MethodDescriptor::ButcherTableauErrorRow[i][MethodDescriptor::StepCount] * tmp.currentEquations[j];
							}
							error[i] += std::pow(sum / errorTemp, 2);
						}
					}
					for (auto& item : error)
					{
						item = std::sqrt(item / problem.DOFCount());
					}

					ValueType stepScale;
					if constexpr (MethodDescriptor::CorrectionMethodsCount == 2)
					{
						tmp.currentError = std::abs(tmp.currentStep) * error[1] / std::sqrt(problem.DOFCount() * (error[1] + 0.01 * error[0]));
						stepScale = this->GetStepScaleFactor() * std::pow(tmp.currentError, -0.7 / MethodDescriptor::AccuracyOrder) * std::pow(tmp.previousError, 0.4 / MethodDescriptor::AccuracyOrder);
					}
					else if constexpr (MethodDescriptor::CorrectionMethodsCount == 1)
					{
						tmp.currentError = error[0];
						stepScale = this->GetStepScaleFactor() * std::pow(tmp.currentError, -1. / (MethodDescriptor::CorrectionMethodsAccuracyOrders[0] + 1));
					}
					if (stepScale > this->GetMaxStepScale())
					{
						stepScale = this->GetMaxStepScale();
					}
					if (stepScale < this->GetMinStepScale())
					{
						stepScale = this->GetMinStepScale();
					}
					tmp.currentStep *= stepScale;

					if (tmp.currentError > this->GetMaxError())
					{
						Logger::Log(MessageType::Info, MessagePriority::Low, MessageTag::TransientSolver,
							"Error overflow on Runge-Kutta step, trying again with lower step...");
						tmp.currentY = std::move(tmp.previousY);
						tmp.currentError = tmp.previousError;
						tmp.currentEquations = std::move(tmp.previousEquations);
						canContinue = false;
					}
				}

				if (canContinue)
				{
					++tmp.stepCount;

					if constexpr (ActualUseAdaptive)
					{
						if (exitConditions & RKExitConditions::StepUnderflow && std::abs(tmp.currentStep) <= std::abs(this->GetMinStep()))
						{
							Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::TransientSolver,
								Format("Stopping Runge-Kutta solver due to step {} decreased below limit {} at time {}.", tmp.currentStep, this->GetMinStep(), tmp.currentX));
							return OutputInfo{ false, OutputInfo::OutputReason::StepUnderflow, tmp.stepCount };
						}
					}
					if (exitConditions & RKExitConditions::SolutionNormOverflow && Norm2(tmp.currentY) >= maxSolutionNorm)
					{
						Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::TransientSolver,
							Format("Stopping Runge-Kutta solver due to solution norm overflow at time {}.", tmp.currentX));
						return OutputInfo{ false, OutputInfo::OutputReason::SolutionNormOverflow, tmp.stepCount };
					}
					if (exitConditions & RKExitConditions::StepCountLimitReached && tmp.stepCount >= stepCountLimit)
					{
						Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::TransientSolver,
							Format("Stopping Runge-Kutta solver due to reaching maximum number of time steps {} at time {}.", stepCountLimit, tmp.currentX));
						return OutputInfo{ false, OutputInfo::OutputReason::StepCountLimitReached, tmp.stepCount };
					}
					tmp.previousX = tmp.currentX;
					tmp.previousEquations = tmp.currentEquations;
					tmp.currentX += tmp.previousStep;

					if (std::abs(lastX - tmp.currentX) < std::abs(tmp.currentStep))
					{
						tmp.currentStep = lastX - tmp.currentX;
					}
				
					if constexpr (ActualUseDenseOutput)
					{
						denseOutput.Process(tmp, this->GetDenseOutputStep(), firstX, problem);
					}
					else
					{
						problem.CacheCurrent();
					}

					Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::TransientSolver,
						Format("Finishing transient problem time step {} at time {} successfully in {}.", tmp.stepCount, tmp.currentX, clock.now() - iterationStartTime));
				}
			}


			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::TransientSolver,
				Format("Finishing transient problem solution successfully after {} steps in {}.", tmp.stepCount, clock.now() - solutionStartTime));

			return OutputInfo{ true, OutputInfo::OutputReason::Success, tmp.stepCount };
		}
	};
}