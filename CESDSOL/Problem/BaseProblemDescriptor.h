#pragma once

#include "Problem/GlobalVariableDependentExpression.h"
#include "Problem/LocalVariableDependentExpression.h"
#include "Problem/Reduction.h"

namespace CESDSOL
{
	template<ProblemType ProblemTypeArg, size_t DimensionArg, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class BaseProblemDescriptor
	{
	public:
		static constexpr size_t Dimension = DimensionArg;
		using CoordinateType = CoordinateTypeArg;
		using FieldType = FieldTypeArg;
		static constexpr ProblemType ProblemType = ProblemTypeArg;
		using GridDescriptorType = GridDescriptor<Dimension, CoordinateType>;

		using CurrentLocalValuesForPIEs = LocalValuesForPIEs<Dimension, ProblemType, CoordinateType, FieldType>;
		using CurrentLocalValuesForVIEs = LocalValuesForVIEs<Dimension, ProblemType, CoordinateType, FieldType>;
		using CurrentLocalValues = LocalValues<Dimension, ProblemType, CoordinateType, FieldType>;
		using CurrentLocalVariableDependentExpression = LocalVariableDependentExpression<Dimension, ProblemType, CoordinateType, FieldType>;
		using CurrentGlobalValues = GlobalValues<ProblemType, CoordinateType, FieldType>;
		using CurrentGlobalValuesForPIEs = GlobalValuesForPIEs<ProblemType, CoordinateType, FieldType>;
		using CurrentGlobalValuesForVIEs = GlobalValuesForVIEs<ProblemType, CoordinateType, FieldType>;
		using CurrentGlobalVariableDependentExpression = GlobalVariableDependentExpression<ProblemType, CoordinateType, FieldType>;
		using CurrentReduction = Reduction<Dimension, ProblemType, CoordinateType, FieldType>;

	protected:
		GridDescriptorType gridDescriptor;
		const Array<Array<std::array<size_t, Dimension>>> derivativeOperators;
		TwoLevelArray<std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>> continuousEquations;
		Array<std::function<FieldType(const CurrentGlobalValues&)>> discreteEquations;
		size_t parameterCount;
		Array<std::function<FieldType(const CurrentLocalValuesForPIEs&,
			const CurrentGlobalValuesForPIEs&)>> localPIEs;
		Array<std::function<FieldType(const CurrentGlobalValuesForPIEs&)>>
			globalPIEs;
		Array<std::function<FieldType(const CurrentLocalValuesForVIEs&,
			const CurrentGlobalValuesForVIEs&)>> localVIEs;
		Array<std::function<FieldType(const CurrentGlobalValuesForVIEs&)>>
			globalVIEs;
		Array<CurrentLocalVariableDependentExpression> localVDEs;
		Array<CurrentGlobalVariableDependentExpression> globalVDEs;
		Array<CurrentReduction> reductions;

		std::string problemName;
		Array<std::string> parameterNames;
		Array<std::string> variableNames;

	public:
		[[nodiscard]] const GridDescriptorType& GetGridDescriptor() const noexcept
		{
			return gridDescriptor;
		}
		
		[[nodiscard]] size_t ParameterCount() const noexcept
		{
			return parameterCount;
		}

		[[nodiscard]] size_t ContinuousEquationCount() const noexcept
		{
			return continuousEquations.size();
		}

		[[nodiscard]] size_t DiscreteEquationCount() const noexcept
		{
			return discreteEquations.size();
		}

		[[nodiscard]] size_t EquationCount() const noexcept
		{
			return ContinuousEquationCount() + DiscreteEquationCount();
		}

		[[nodiscard]] size_t DerivativeOperatorCount(size_t fieldIndex) const noexcept
		{
			return derivativeOperators[fieldIndex].size();
		}

		[[nodiscard]] const std::array<size_t, Dimension>& GetDerivativeOperator(size_t fieldIndex, size_t operatorIndex) const noexcept
		{
			return derivativeOperators[fieldIndex][operatorIndex];
		}

		[[nodiscard]] size_t LocalPIECount() const noexcept
		{
			return localPIEs.size();
		}

		[[nodiscard]] size_t GlobalPIECount() const noexcept
		{
			return globalPIEs.size();
		}

		[[nodiscard]] size_t LocalVIECount() const noexcept
		{
			return localVIEs.size();
		}

		[[nodiscard]] size_t GlobalVIECount() const noexcept
		{
			return globalVIEs.size();
		}

		[[nodiscard]] size_t LocalVDECount() const noexcept
		{
			return localVDEs.size();
		}

		[[nodiscard]] size_t GlobalVDECount() const noexcept
		{
			return globalVDEs.size();
		}

		[[nodiscard]] size_t ReductionCount() const noexcept
		{
			return reductions.size();
		}

		void SetContinuousEquation(
			size_t equationIndex,
			size_t regionIndex,
			std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> equation
		) noexcept
		{
			continuousEquations[equationIndex][regionIndex] = equation;
		}

		void SetDiscreteEquation(
			size_t equationIndex,
			std::function<FieldType(const CurrentGlobalValues&)> equation
		) noexcept
		{
			discreteEquations[equationIndex] = equation;
		}

		void SetLocalParameterIndependentExpression(
			size_t expressionIndex,
			std::function<FieldType(const CurrentLocalValuesForPIEs&,
				const CurrentGlobalValuesForPIEs&)> expression
		) noexcept
		{
			localPIEs[expressionIndex] = expression;
		}

		void SetGlobalParameterIndependentExpression(
			size_t expressionIndex,
			std::function<FieldType(const CurrentGlobalValuesForPIEs&)> expression
		) noexcept
		{
			globalPIEs[expressionIndex] = expression;
		}

		void SetLocalVariableIndependentExpression(
			size_t expressionIndex,
			std::function<FieldType(const CurrentLocalValuesForVIEs&,
				const CurrentGlobalValuesForVIEs&)> expression
		) noexcept
		{
			localVIEs[expressionIndex] = expression;
		}

		void SetGlobalVariableIndependentExpression(
			size_t expressionIndex,
			std::function<FieldType(const CurrentGlobalValuesForVIEs&)> expression
		) noexcept
		{
			globalVIEs[expressionIndex] = expression;
		}

		void SetLocalVariableDependentExpression(
			size_t expressionIndex,
			const std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>& expression
		) noexcept
		{
			localVDEs[expressionIndex].Expression = expression;
		}

		void SetGlobalVariableDependentExpression(
			size_t expressionIndex,
			const std::function<FieldType(const CurrentGlobalValues&)>& expression
		) noexcept
		{
			globalVDEs[expressionIndex].Expression = expression;
		}

		void SetReduction(
			size_t reductionIndex,
			const std::function<FieldType(FieldType)>& externalFunction,
			const std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>& internalFunction
		) noexcept
		{
			reductions[reductionIndex].ExternalFunction = externalFunction;
			reductions[reductionIndex].InternalFunction = internalFunction;
		}

		void SetIntegrand(
			size_t reductionIndex,
			const std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>& integrand
		) noexcept
		{
			reductions[reductionIndex].InternalFunction =
				[=](const auto& locals, const auto& globals) {return locals.IntegrationWeight * integrand(locals, globals); };
		}

		[[nodiscard]] FieldType CalculateContinuousEquation(
			size_t equationIndex,
			size_t regionIndex,
			const CurrentLocalValues& locals,
			const CurrentGlobalValues& globals
		) const noexcept
		{
			return continuousEquations[equationIndex][regionIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateDiscreteEquation(
			size_t equationIndex,
			const CurrentGlobalValues& globals
		) const noexcept
		{
			return discreteEquations[equationIndex](globals);
		}

		[[nodiscard]] FieldType CalculateLocalParameterIndependentExpression(
			size_t expressionIndex,
			const CurrentLocalValuesForPIEs& locals,
			const CurrentGlobalValuesForPIEs& globals
		) const noexcept
		{
			return localPIEs[expressionIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateGlobalParameterIndependentExpression(
			size_t expressionIndex,
			const CurrentGlobalValuesForPIEs& globals
		) const noexcept
		{
			return globalPIEs[expressionIndex](globals);
		}

		[[nodiscard]] FieldType CalculateLocalVariableIndependentExpression(
			size_t expressionIndex,
			const CurrentLocalValuesForVIEs& locals,
			const CurrentGlobalValuesForVIEs& globals
		) const noexcept
		{
			return localVIEs[expressionIndex](locals, globals);
		}

		[[nodiscard]] FieldType CalculateGlobalVariableIndependentExpression(
			size_t expressionIndex,
			const CurrentGlobalValuesForVIEs& globals
		) const noexcept
		{
			return globalVIEs[expressionIndex](globals);
		}

		[[nodiscard]] FieldType CalculateLocalVariableDependentExpression(
			size_t expressionIndex,
			const CurrentLocalValues& locals,
			const CurrentGlobalValues& globals
		) const noexcept
		{
			return localVDEs[expressionIndex].Expression(locals, globals);
		}

		[[nodiscard]] FieldType CalculateGlobalVariableDependentExpression(
			size_t expressionIndex,
			const CurrentGlobalValues& globals) const noexcept
		{
			return globalVDEs[expressionIndex].Expression(globals);
		}

		[[nodiscard]] FieldType CalculateReductionPoint(
			size_t reductionIndex,
			const CurrentLocalValues& locals,
			const CurrentGlobalValues& globals
		) const noexcept
		{
			return reductions[reductionIndex].InternalFunction(locals, globals);
		}

		[[nodiscard]] FieldType CalculateReductionTotal(
			size_t reductionIndex,
			FieldType sum
		) const noexcept
		{
			return reductions[reductionIndex].ExternalFunction(sum);
		}

		[[nodiscard]] bool HasContinuousEquation(size_t equationIndex, size_t regionIndex) const noexcept
		{
			return static_cast<bool>(continuousEquations[equationIndex][regionIndex]);
		}

		void SetProblemName(const std::string& name) noexcept
		{
			problemName = name;
		}

		[[nodiscard]] const std::string& GetProblemName() const noexcept
		{
			return problemName;
		}

		void SetParameterName(size_t index, const std::string& name) noexcept
		{
			parameterNames.At(index) = name;
		}

		void SetVariableName(size_t index, const std::string& name) noexcept
		{
			variableNames.At(index) = name;
		}

		[[nodiscard]] const std::string& GetParameterName(size_t index) const noexcept
		{
			return parameterNames.At(index);
		}

		[[nodiscard]] const std::string& GetVariableName(size_t index) const noexcept
		{
			return variableNames.At(index);
		}

		BaseProblemDescriptor(
			const GridDescriptorType& aGridDescriptor,
			const Array<Array<std::array<size_t, Dimension>>>& aDerivativeOperators,
			size_t continuousEquationsCount = 1,
			size_t aParameterCount = 0,
			size_t discreteEquationsCount = 0,
			size_t localParameterIndependentExpressionsCount = 0,
			size_t globalParameterIndependentExpressionsCount = 0,
			size_t localVariableIndependentExpressionsCount = 0,
			size_t globalVariableIndependentExpressionsCount = 0,
			size_t localVariableDependentExpressionsCount = 0,
			size_t globalVariableDependentExpressionsCount = 0,
			size_t reductionsCount = 0)
			: gridDescriptor(aGridDescriptor)
			, derivativeOperators(aDerivativeOperators)
			, continuousEquations({ {continuousEquationsCount, gridDescriptor.GetRegionCount()} })
			, discreteEquations(discreteEquationsCount)
			, parameterCount(aParameterCount)
			, localPIEs(localParameterIndependentExpressionsCount)
			, globalPIEs(globalParameterIndependentExpressionsCount)
			, localVIEs(localVariableIndependentExpressionsCount)
			, globalVIEs(globalVariableIndependentExpressionsCount)
			, localVDEs(localVariableDependentExpressionsCount)
			, globalVDEs(globalVariableDependentExpressionsCount)
			, reductions(reductionsCount)
			, parameterNames(parameterCount)
			, variableNames(EquationCount())
		{}

		BaseProblemDescriptor(
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
			: BaseProblemDescriptor(
				aGridDescriptor,
				Array(derivativeOperators, continuousEquationsCount),
				continuousEquationsCount,
				parameterCount,
				discreteEquationsCount,
				localParameterIndependentExpressionsCount,
				globalParameterIndependentExpressionsCount,
				localVariableIndependentExpressionsCount,
				globalVariableIndependentExpressionsCount,
				localVariableDependentExpressionsCount,
				globalVariableDependentExpressionsCount,
				reductionsCount)
		{}

		[[nodiscard]] virtual bool Validate() const noexcept
		{
			const auto andAccumulator = [](auto sum, const auto& element) {return sum && element; };
			const auto validityAndAccumulator = [](auto sum, const auto& element) {return sum && element.IsValid(); };
			return std::accumulate(continuousEquations.begin(), continuousEquations.end(), true, [](auto sum, const auto& element) {return sum && element[0]; })
				|| std::accumulate(discreteEquations.begin(), discreteEquations.end(), true, andAccumulator)
				|| std::accumulate(localPIEs.begin(), localPIEs.end(), true, andAccumulator)
				|| std::accumulate(globalPIEs.begin(), globalPIEs.end(), true, andAccumulator)
				|| std::accumulate(localVIEs.begin(), localVIEs.end(), true, andAccumulator)
				|| std::accumulate(globalVIEs.begin(), globalVIEs.end(), true, andAccumulator)
				|| std::accumulate(localVDEs.begin(), localVDEs.end(), true, validityAndAccumulator)
				|| std::accumulate(globalVDEs.begin(), globalVDEs.end(), true, validityAndAccumulator)
				|| std::accumulate(reductions.begin(), reductions.end(), true, validityAndAccumulator);
		}

		~BaseProblemDescriptor() = default;
	};
}