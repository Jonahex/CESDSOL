#pragma once

#include "Discretization/Discretization.h"
#include "Grid/Grid.h"
#include "Problem/BaseProblemDescriptor.h"

namespace CESDSOL
{
	template<size_t DimensionArg, template<typename> typename MatrixTypeArg = CSRMatrix, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class StationaryProblem;
	
	template<size_t DimensionArg, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class StationaryProblemDescriptor
		: public BaseProblemDescriptor<ProblemType::Stationary, DimensionArg, CoordinateTypeArg, FieldTypeArg>
	{
	private:
		using BaseType = BaseProblemDescriptor<ProblemType::Stationary, DimensionArg, CoordinateTypeArg, FieldTypeArg>;
		
	public:
		using typename BaseType::CoordinateType;
		using typename BaseType::FieldType;
		using typename BaseType::GridDescriptorType;

		using BaseType::DiscreteEquationCount;
		using BaseType::ContinuousEquationCount;
		using BaseType::DerivativeOperatorCount;
		using BaseType::LocalVDECount;
		using BaseType::ReductionCount;

		using BaseType::gridDescriptor;
		using BaseType::continuousEquations;
		using BaseType::discreteEquations;
		using BaseType::localPIEs;
		using BaseType::globalPIEs;
		using BaseType::localVIEs;
		using BaseType::globalVIEs;
		using BaseType::globalVDEs;
		using BaseType::localVDEs;
		using BaseType::reductions;
		
		static constexpr size_t Dimension = BaseType::Dimension;
		static constexpr ProblemType ProblemType = BaseType::ProblemType;
		
		using CurrentLocalValuesForJacobian = LocalValuesForJacobian<Dimension, ProblemType, CoordinateType, FieldType>;
		using CurrentGlobalValuesForJacobian = GlobalValuesForJacobian<ProblemType, CoordinateType, FieldType>;

	private:
		[[nodiscard]] static FieldType DefaultMeritFunction(const Array<FieldType>& fields) noexcept
		{
			return Norm2(fields) / fields.size();
		}

		FourLevelArray<std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>> jacobianComponents;
		std::function<FieldType(const Array<FieldType>&)> meritFunction = DefaultMeritFunction;

		[[nodiscard]] auto ConstructContinuousEquationJacobianLevelStructure() const noexcept
		{
			const bool hasDiscrete = DiscreteEquationCount() > 0;
			Array<std::pair<size_t, LevelStructure<2>>> result{ ContinuousEquationCount() + static_cast<bool>(hasDiscrete) };
			for (size_t i = 0; i < ContinuousEquationCount(); i++)
			{
				result[i] = { 1, {DerivativeOperatorCount(i) + 1, gridDescriptor.GetRegionCount()} };
			}
			if (hasDiscrete)
			{
				result[ContinuousEquationCount()] = { DiscreteEquationCount(), {1, gridDescriptor.GetRegionCount()} };
			}
			return result;
		}

		[[nodiscard]] auto ConstructDiscreteEquationJacobianLevelStructure() const noexcept
		{
			Array<std::pair<size_t, LevelStructure<2>>> result{ ContinuousEquationCount() + 1 };
			for (size_t i = 0; i < ContinuousEquationCount(); i++)
			{
				result[i] = { 1, {DerivativeOperatorCount(i) + 1, 1} };
			}
			result[ContinuousEquationCount()] = { DiscreteEquationCount(), {1, 1} };
			return result;
		}

		[[nodiscard]] auto ConstructJacobianComponents() const noexcept
		{
			if (DiscreteEquationCount() > 0)
			{
				return FourLevelArray<std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>>
				({ {ContinuousEquationCount(), ConstructContinuousEquationJacobianLevelStructure()}, {DiscreteEquationCount(), ConstructDiscreteEquationJacobianLevelStructure()} });
			}
			else
			{
				return FourLevelArray<std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>>
					({ContinuousEquationCount(), ConstructContinuousEquationJacobianLevelStructure() });
			}
		}

		void InitNestedJacobians() noexcept
		{
			for (auto& gvde : globalVDEs)
			{
				gvde.JacobianComponents = { DiscreteEquationCount() };
			}
			
			if (LocalVDECount() == 0 && ReductionCount() == 0)
			{
				return;
			}
			
			const bool hasDiscrete = DiscreteEquationCount() > 0;
			Array<std::pair<size_t, size_t>> levelStructure(ContinuousEquationCount() + static_cast<bool>(hasDiscrete));
			for (size_t i = 0; i < ContinuousEquationCount(); i++)
			{
				levelStructure[i] = { 1, DerivativeOperatorCount(i) + 1 };
			}
			if (hasDiscrete)
			{
				levelStructure[ContinuousEquationCount()] = { DiscreteEquationCount(), 1 };
			}
			
			for (auto& lvde : localVDEs)
			{
				lvde.JacobianComponents = {levelStructure};
			}
			for (auto& reduction : reductions)
			{
				reduction.InternalFunctionJacobianComponents = { levelStructure };
			}
		}

	public:
		void SetJacobianComponent(
			size_t equationIndex, 
			size_t fieldIndex, 
			size_t operatorIndex, 
			size_t regionIndex, 
			std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)> component
		) noexcept
		{
			jacobianComponents[equationIndex][fieldIndex][operatorIndex][regionIndex] = component;
		}
		
		void SetLocalVariableDependentExpressionJacobianComponent(
			size_t expressionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>& expression
		) noexcept
		{
			localVDEs[expressionIndex].JacobianComponents[fieldIndex][operatorIndex] = expression;
		}

		void SetGlobalVariableDependentExpressionJacobianComponent(
			size_t expressionIndex,
			size_t variableIndex,
			const std::function<FieldType(const CurrentGlobalValuesForJacobian&)>& expression
		) noexcept
		{
			globalVDEs[expressionIndex].JacobianComponents[variableIndex] = expression;
		}

		void SetReductionExternalJacobian(
			size_t reductionIndex,
			const std::function<FieldType(FieldType)>& expression
		) noexcept
		{
			reductions[reductionIndex].ExternalFunctionJacobian = expression;
		}

		void SetReductionInternalJacobianComponent(
			size_t reductionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>& expression
		) noexcept
		{
			reductions[reductionIndex].InternalFunctionJacobianComponents[fieldIndex][operatorIndex] = expression;
		}

		void SetIntegrandJacobianComponent(
			size_t reductionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const std::function<FieldType(const CurrentLocalValuesForJacobian&, const CurrentGlobalValuesForJacobian&)>& integrand
		) noexcept
		{
			reductions[reductionIndex].InternalFunctionJacobianComponents[fieldIndex][operatorIndex] = 
				[=](const auto& locals, const auto& globals) {return locals.IntegrationWeight * integrand(locals, globals); };
		}

		void SetMeritFunction(std::function<FieldType(const Array<FieldType>&)> function) noexcept
		{
			meritFunction = function;
		}

		[[nodiscard]] FieldType CalculateJacobianComponent(
			size_t equationIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			size_t regionIndex,
			const CurrentLocalValuesForJacobian& locals,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return jacobianComponents[equationIndex][fieldIndex][operatorIndex][regionIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateLVDEJacobianComponent(
			size_t expressionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const CurrentLocalValuesForJacobian& locals,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return localVDEs[expressionIndex].JacobianComponents[fieldIndex][operatorIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateGVDEJacobianComponent(
			size_t expressionIndex,
			size_t fieldIndex,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return globalVDEs[expressionIndex].JacobianComponents[fieldIndex](globals);
		}

		[[nodiscard]] FieldType CalculateReductionExternalJacobianComponent(
			size_t reductionIndex,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return reductions[reductionIndex].ExternalFunctionJacobian(globals.ReductionValues[reductionIndex]);
		}

		[[nodiscard]] FieldType CalculateReductionInternalJacobianComponent(
			size_t reductionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const CurrentLocalValuesForJacobian& locals,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return reductions[reductionIndex].InternalFunctionJacobianComponents[fieldIndex][operatorIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateReductionJacobianComponent(
			size_t reductionIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			const CurrentLocalValuesForJacobian& locals,
			const CurrentGlobalValuesForJacobian& globals
		) const noexcept
		{
			return reductions[reductionIndex].ExternalFunctionJacobian(globals.ReductionValues[reductionIndex]) *
				reductions[reductionIndex].InternalFunctionJacobianComponents[fieldIndex][operatorIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateMerit(const Array<FieldType>& variables) const noexcept
		{
			return meritFunction(variables);
		}

		[[nodiscard]] bool HasJacobianComponent(
			size_t equationIndex,
			size_t fieldIndex,
			size_t operatorIndex,
			size_t regionIndex
		) const noexcept
		{
			return static_cast<bool>(jacobianComponents[equationIndex][fieldIndex][operatorIndex][regionIndex]);
		}

		[[nodiscard]] bool HasReductionJacobianComponent(
			size_t reductionIndex,
			size_t fieldIndex,
			size_t operatorIndex
		) const noexcept
		{
			return static_cast<bool>(reductions[reductionIndex].InternalFunctionJacobianComponents[fieldIndex][operatorIndex]);
		}

		[[nodiscard]] bool HasLVDEJacobianComponent(
			size_t expressionIndex,
			size_t fieldIndex,
			size_t operatorIndex
		) const noexcept
		{
			return static_cast<bool>(localVDEs[expressionIndex].JacobianComponents[fieldIndex][operatorIndex]);
		}

		[[nodiscard]] bool HasGVDEJacobianComponent(
			size_t expressionIndex,
			size_t variableIndex
		) const noexcept
		{
			return static_cast<bool>(globalVDEs[expressionIndex].JacobianComponents[variableIndex]);
		}

		template<template<typename> typename MatrixType = CSRMatrix>
		[[nodiscard]] auto
			MakeProblem(sptr<Grid<Dimension, CoordinateType>> grid, uptr<Discretization<Dimension, MatrixType, CoordinateType>> discretizer) const noexcept;

		StationaryProblemDescriptor(
			const GridDescriptorType& aGridDescriptor,
			const Array<Array<std::array<size_t, Dimension>>>& derivativeOperators,
			size_t continuousEquationsCount = 1,
			size_t parameterCount = 0,
			size_t discreteEquationsCount = 0,
			size_t localParameterIndependentExpressionsCount = 0,
			size_t globalParameterIndependentExpressionsCount = 0,
			size_t localVariableIndependentExpressionsCount = 0,
			size_t globalVariableIndependentExpressionsCount = 0,
			size_t localVariableDependentExpressionsCount = 0,
			size_t globalVariableDependentExpressionsCount = 0,
			size_t reductionsCount = 0)
			: BaseType
			( aGridDescriptor
			, derivativeOperators
			, continuousEquationsCount
			, parameterCount
			, discreteEquationsCount
			, localParameterIndependentExpressionsCount
			, globalParameterIndependentExpressionsCount
			, localVariableIndependentExpressionsCount
			, globalVariableIndependentExpressionsCount
			, localVariableDependentExpressionsCount
			, globalVariableDependentExpressionsCount
			, reductionsCount
			)
			, jacobianComponents(ConstructJacobianComponents())
		{
			InitNestedJacobians();
		}

		StationaryProblemDescriptor(
			const GridDescriptorType& aGridDescriptor,
			const Array<std::array<size_t, Dimension>>& derivativeOperators,
			size_t continuousEquationsCount = 1,
			size_t parameterCount = 0,
			size_t discreteEquationsCount = 0,
			size_t localParameterIndependentExpressionsCount = 0,
			size_t globalParameterIndependentExpressionsCount = 0,
			size_t localVariableIndependentExpressionsCount = 0,
			size_t globalVariableIndependentExpressionsCount = 0,
			size_t localVariableDependentExpressionsCount = 0,
			size_t globalVariableDependentExpressionsCount = 0,
			size_t reductionsCount = 0)
			: StationaryProblemDescriptor
			( aGridDescriptor
			, Array(derivativeOperators, continuousEquationsCount)
			, continuousEquationsCount
			, parameterCount
			, discreteEquationsCount
			, localParameterIndependentExpressionsCount
			, globalParameterIndependentExpressionsCount
			, localVariableIndependentExpressionsCount
			, globalVariableIndependentExpressionsCount
			, localVariableDependentExpressionsCount
			, globalVariableDependentExpressionsCount
			, reductionsCount
			)
		{}
	};
}