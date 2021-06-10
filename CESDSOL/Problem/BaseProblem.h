#pragma once

#include "Discretization/Discretization.h"
#include "Grid/Grid.h"
#include "Problem/BaseProblemDescriptor.h"
#include "Problem/ProblemUtils.h"
#include "Serialization/DataToLoad.h"
#include "Serialization/DataType.h"
#include "Serialization/ProblemType.h"

#include <any>

namespace CESDSOL
{
	template<typename DescriptorTypeArg, template<typename> typename MatrixTypeArg = CSRMatrix>
	class BaseProblem
	{
	private:
		using Self = BaseProblem<DescriptorTypeArg, MatrixTypeArg>;

	public:
		using DescriptorType = DescriptorTypeArg;
		static constexpr size_t Dimension = DescriptorType::Dimension;
		using CoordinateType = DescriptorType::CoordinateType;
		using DiscretizationType = Discretization<Dimension, MatrixTypeArg, CoordinateType>;
		using FieldType = DescriptorType::FieldType;
		using GridType = Grid<Dimension, CoordinateType>;
		static constexpr ProblemType ProblemType = DescriptorType::ProblemType;

		using CurrentLocalValuesForPIEs = typename DescriptorType::CurrentLocalValuesForPIEs;
		using CurrentLocalValuesForVIEs = typename DescriptorType::CurrentLocalValuesForVIEs;
		using CurrentLocalValues = typename DescriptorType::CurrentLocalValues;
		using CurrentGlobalValuesForPIEs = typename DescriptorType::CurrentGlobalValuesForPIEs;
		using CurrentGlobalValuesForVIEs = typename DescriptorType::CurrentGlobalValuesForVIEs;
		using CurrentGlobalValues = typename DescriptorType::CurrentGlobalValues;
		using CurrentReduction = typename DescriptorType::CurrentReduction;

		using DifferentiationMatrixType = MatrixTypeArg<CoordinateType>;
		using VectorType = Vector<FieldType>;

	protected:
		sptr<GridType> grid;
		uptr<DiscretizationType> discretizer;
		DescriptorType descriptor;

		Array<std::array<size_t, Dimension>> derivativeOperators;
		TwoLevelArray<size_t> fieldDerivativeOperatorMap;

		Array<FieldType> parameters;
		TwoLevelArray<FieldType> variables;
		ThreeLevelArray<FieldType> derivatives;
		TwoLevelArray<FieldType> equations;

		TwoLevelArray<FieldType> localPIEs;
		Array<FieldType> globalPIEs;
		TwoLevelArray<FieldType> localVIEs;
		Array<FieldType> globalVIEs;
		TwoLevelArray<FieldType> localVDEs;
		Array<FieldType> globalVDEs;
		Array<FieldType> reductions;

		Array<DifferentiationMatrixType> differentiationWeights;
		Array<CoordinateType> integrationWeights;

		std::vector<std::pair<std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>, std::string>> localOutputExpressions;
		std::vector<std::pair<std::function<FieldType(const CurrentGlobalValues&)>, std::string>> globalOutputExpressions;
		std::vector<std::pair<CurrentReduction, std::string>> reductionOutputExpressions;

		struct PointOutputExpression
		{
			std::string name;
			std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> expression;
			std::array<CoordinateType, Dimension> point;
			SparseVector<CoordinateType> weights;
		};

		std::vector<PointOutputExpression> pointOutputExpressions;

		bool isActualOnParameters = true;
		bool isActualOnVariables = true;

	public:
		[[nodiscard]] size_t DOFCount() const noexcept
		{
			return grid->GetSize() * descriptor.ContinuousEquationCount() + descriptor.DiscreteEquationCount();
		}

	protected:
		[[nodiscard]] void EnumerateDerivativeOperators() noexcept
		{
			std::unordered_map<std::array<size_t, Dimension>, size_t> operators;
			size_t count = 0;
			Array<std::pair<size_t, size_t>> levelStructure{ descriptor.ContinuousEquationCount() };
			for (size_t i = 0; i < descriptor.ContinuousEquationCount(); i++)
			{
				levelStructure[i] = { 1, descriptor.DerivativeOperatorCount(i) };
				for (size_t j = 0; j < descriptor.DerivativeOperatorCount(i); j++)
				{
					const auto& currentOperator = descriptor.GetDerivativeOperator(i, j);
					if (!operators.contains(currentOperator))
					{
						operators.insert({ currentOperator, count++ });
					}
				}
			}
			fieldDerivativeOperatorMap = TwoLevelArray<size_t>(levelStructure);
			derivativeOperators = Array<std::array<size_t, Dimension>>(operators.size());
			for (const auto& x : operators)
			{
				derivativeOperators[x.second] = x.first;
			}
			for (size_t i = 0; i < descriptor.ContinuousEquationCount(); i++)
			{
				for (size_t j = 0; j < descriptor.DerivativeOperatorCount(i); j++)
				{
					fieldDerivativeOperatorMap[i][j] = operators[descriptor.GetDerivativeOperator(i, j)];
				}
			}
		}

		[[nodiscard]] ThreeLevelArray<FieldType> ConstructDerivatives() const noexcept
		{
			Array<std::pair<size_t, LevelStructure<2>>> levelStructure{ descriptor.ContinuousEquationCount() };
			for (size_t i = 0; i < descriptor.ContinuousEquationCount(); i++)
			{
				levelStructure[i] = { 1, { descriptor.DerivativeOperatorCount(i), grid->GetSize() } };
			}
			return ThreeLevelArray<FieldType>(levelStructure);
		}

		void ConstructDifferentiationWeights() noexcept
		{
			differentiationWeights = Array<CSRMatrix<CoordinateType>>(derivativeOperators.size());
			for (size_t i = 0; i < derivativeOperators.size(); i++)
			{
				differentiationWeights[i] = discretizer->GetDifferentiationMatrix(*grid, derivativeOperators[i]);
			}
		}

		void ConstructIntegrationWeights() noexcept
		{
			integrationWeights = discretizer->GetIntegrationWeightsVector(*grid);
		}

		void CalculateParameterIndependentExpressions() noexcept
		{
			CurrentLocalValuesForPIEs locals{ descriptor.LocalPIECount() };
			CurrentGlobalValuesForPIEs globals{ globalPIEs };
			for (size_t i = 0; i < descriptor.GlobalPIECount(); i++)
			{
				globalPIEs[i] = descriptor.CalculateGlobalParameterIndependentExpression(i, globals);
			}
			for (size_t i = 0; i < grid->GetSize(); i++)
			{
				locals.Point = grid->GetCoordinates(i);
				for (size_t j = 0; j < descriptor.LocalPIECount(); j++)
				{
					localPIEs[j][i] = descriptor.CalculateLocalParameterIndependentExpression(j, locals, globals);
					locals.PIEValues[i] = localPIEs[j][i];
				}
			}
		}

		[[nodiscard]] Array<std::pair<size_t, size_t>> ConstructDerivativesLevelStructure() const noexcept
		{
			Array<std::pair<size_t, size_t>> derivativesLevelStructure{ descriptor.ContinuousEquationCount() };
			for (size_t i = 0; i < derivativesLevelStructure.size(); i++)
			{
				derivativesLevelStructure[i] = { 1, descriptor.DerivativeOperatorCount(i) };
			}
			return derivativesLevelStructure;
		}

		CurrentLocalValuesForVIEs ConstructLocalValuesForVIEs() const noexcept
		{
			return CurrentLocalValuesForVIEs
			{
				descriptor.LocalPIECount(),
				descriptor.LocalVIECount()
			};
		}

		CurrentLocalValues ConstructLocalValues() const noexcept
		{
			return CurrentLocalValues
			{
				descriptor.LocalPIECount(),
				descriptor.LocalVIECount(),
				descriptor.ContinuousEquationCount(),
				ConstructDerivativesLevelStructure(),
				descriptor.LocalVDECount()
			};
		}

		virtual CurrentGlobalValuesForVIEs ConstructGlobalValuesForVIEs() const noexcept
		{
			return CurrentGlobalValuesForVIEs
			{
				globalPIEs,
				parameters,
				globalVIEs
			};
		}

		CurrentGlobalValues ConstructGlobalValues() const noexcept
		{
			return CurrentGlobalValues
			{
				globalPIEs,
				parameters,
				globalVIEs,
				std::span(variables[descriptor.ContinuousEquationCount()].data(), variables[descriptor.ContinuousEquationCount()].data() + descriptor.DiscreteEquationCount()),
				globalVDEs,
				reductions
			};
		}

		void FillEssentialLocals(size_t pointIndex, CurrentLocalValues& locals) const noexcept
		{
			locals.Point = grid->GetCoordinates(pointIndex);
			locals.IntegrationWeight = integrationWeights[pointIndex];
			for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
			{
				locals.FieldValues[j] = variables[j][pointIndex];
				for (size_t k = 0; k < descriptor.DerivativeOperatorCount(j); k++)
				{
					locals.DerivativeValues[j][k] = derivatives[j][k][pointIndex];
				}
			}
		}

		void FillPIEs(size_t pointIndex, CurrentLocalValuesForPIEs& locals) const noexcept
		{
			for (size_t j = 0; j < descriptor.LocalPIECount(); j++)
			{
				locals.PIEValues[j] = localPIEs[j][pointIndex];
			}
		}

		void FillVIEs(size_t pointIndex, CurrentLocalValuesForVIEs& locals) const noexcept
		{
			for (size_t j = 0; j < descriptor.LocalVIECount(); j++)
			{
				locals.VIEValues[j] = localVIEs[j][pointIndex];
			}
		}

		void FillVDEs(size_t pointIndex, CurrentLocalValues& locals) const noexcept
		{
			for (size_t j = 0; j < descriptor.LocalVDECount(); j++)
			{
				locals.VDEValues[j] = localVDEs[j][pointIndex];
			}
		}

		void FillAllLocals(size_t pointIndex, CurrentLocalValues& locals) const noexcept
		{
			FillEssentialLocals(pointIndex, locals);
			FillPIEs(pointIndex, locals);
			FillVIEs(pointIndex, locals);
			FillVDEs(pointIndex, locals);
		}

		[[nodiscard]] size_t GetTrueRegionIndex(size_t equationIndex, size_t pointIndex) const noexcept
		{
			const auto regionIndex = grid->GetRegionIndex(pointIndex);
			return descriptor.HasContinuousEquation(equationIndex, regionIndex) ? regionIndex : 0;
		}

		void UpdateDerivatives() noexcept
		{
			for (size_t i = 0; i < descriptor.ContinuousEquationCount(); i++)
			{
				for (size_t j = 0; j < descriptor.DerivativeOperatorCount(i); j++)
				{
					Multiply(differentiationWeights[fieldDerivativeOperatorMap[i][j]], variables[i], derivatives[i][j]);
				}
			}
		}

		void UpdateVariableIndependentExpressions() noexcept
		{
			auto globals = ConstructGlobalValuesForVIEs();
			for (size_t i = 0; i < descriptor.GlobalVIECount(); i++)
			{
				globalVIEs[i] = descriptor.CalculateGlobalVariableIndependentExpression(i, globals);
			}
#pragma omp parallel
			{
				auto locals = ConstructLocalValuesForVIEs();
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); i++)
				{
					locals.Point = grid->GetCoordinates(i);
					FillPIEs(i, locals);
					for (size_t j = 0; j < descriptor.LocalVIECount(); j++)
					{
						localVIEs[j][i] = descriptor.CalculateLocalVariableIndependentExpression(j, locals, globals);
						locals.VIEValues[j] = localVIEs[j][i];
					}
				}
			}
		}

		void UpdateReductions() noexcept
		{
			auto locals = ConstructLocalValues();
			const auto globals = ConstructGlobalValues();
			std::fill(reductions.begin(), reductions.end(), 0);
			for (size_t i = 0; i < grid->GetSize(); i++)
			{
				FillAllLocals(i, locals);
				for (size_t j = 0; j < descriptor.ReductionCount(); j++)
				{
					reductions[j] += descriptor.CalculateReductionPoint(j, locals, globals);
				}
			}
			for (size_t j = 0; j < descriptor.ReductionCount(); j++)
			{
				reductions[j] = descriptor.CalculateReductionTotal(j, reductions[j]);
			}
		}

		void UpdateVariableDependentExpressions() noexcept
		{
			auto globals = ConstructGlobalValues();
			for (size_t i = 0; i < descriptor.GlobalVDECount(); i++)
			{
				globalVDEs[i] = descriptor.CalculateGlobalVariableDependentExpression(i, globals);
			}
#pragma omp parallel
			{
				auto locals = ConstructLocalValues();
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); i++)
				{
					FillEssentialLocals(i, locals);
					FillPIEs(i, locals);
					FillVIEs(i, locals);
					for (size_t j = 0; j < descriptor.LocalVDECount(); j++)
					{
						localVDEs[j][i] = descriptor.CalculateLocalVariableDependentExpression(j, locals, globals);
						locals.VDEValues[j] = localVDEs[j][i];
					}
				}
			}
		}

		void UpdateEquations() noexcept
		{
			const auto globals = ConstructGlobalValues();
			for (size_t i = 0; i < descriptor.DiscreteEquationCount(); i++)
			{
				equations[descriptor.ContinuousEquationCount() + i][0] = descriptor.CalculateDiscreteEquation(i, globals);
			}
#pragma omp parallel
			{
				auto locals = ConstructLocalValues();
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); i++)
				{
					FillAllLocals(i, locals);
					const auto regionIndex = grid->GetRegionIndex(i);
					for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
					{
						const auto trueRegionIndex = descriptor.HasContinuousEquation(j, regionIndex) ? regionIndex : 0;
						equations[j][i] = descriptor.CalculateContinuousEquation(j, trueRegionIndex, locals, globals);
					}
				}
			}
		}

		void Actualize() noexcept
		{
			if (!isActualOnVariables)
			{
				UpdateDerivatives();
			}
			if (!isActualOnParameters)
			{
				UpdateVariableIndependentExpressions();
			}
			if (!isActualOnVariables || !isActualOnParameters)
			{
				UpdateVariableDependentExpressions();
				UpdateReductions();
				UpdateEquations();
			}
			isActualOnVariables = true;
			isActualOnParameters = true;
		}

		BaseProblem(sptr<GridType> aGrid, uptr<DiscretizationType> aDiscretizer, const DescriptorType& descriptor)
			: grid(std::move(aGrid))
			, discretizer(std::move(aDiscretizer))
			, descriptor(descriptor)
			, parameters(descriptor.ParameterCount())
			, variables({ {descriptor.ContinuousEquationCount(), grid->GetSize()}, {descriptor.DiscreteEquationCount(), 1} })
			, derivatives(ConstructDerivatives())
			, equations({ {descriptor.ContinuousEquationCount(), grid->GetSize()}, {descriptor.DiscreteEquationCount(), 1} })
			, localPIEs({ descriptor.LocalPIECount(), grid->GetSize() })
			, globalPIEs(descriptor.GlobalPIECount())
			, localVIEs({ descriptor.LocalVIECount(), grid->GetSize() })
			, globalVIEs(descriptor.GlobalVIECount())
			, localVDEs({ descriptor.LocalVDECount(), grid->GetSize() })
			, globalVDEs(descriptor.GlobalVDECount())
			, reductions(descriptor.ReductionCount())
		{
			EnumerateDerivativeOperators();
			ConstructDifferentiationWeights();
			ConstructIntegrationWeights();
			CalculateParameterIndependentExpressions();
		}

	public:
		[[nodiscard]] const GridType& GetGrid() const noexcept
		{
			return *grid;
		}

		[[nodiscard]] const DiscretizationType& GetDiscretizer() const noexcept
		{
			return *discretizer;
		}

		[[nodiscard]] const DescriptorType& GetDescriptor() const noexcept
		{
			return descriptor;
		}

		[[nodiscard]] size_t ParameterCount() const noexcept
		{
			return descriptor.ParameterCount();
		}

		void SetParameter(size_t parameterIndex, FieldType value) noexcept
		{
			parameters[parameterIndex] = value;
			isActualOnParameters = false;
		}

		void SetParameters(const Vector<FieldType>& aParameters) noexcept
		{
			AssertE(aParameters.size() == ParameterCount(), MessageTag::Problem, "Trying to set inconsistent number of parameters!");
			parameters = aParameters;
			isActualOnParameters = false;
		}

		[[nodiscard]] FieldType GetParameter(size_t parameterIndex) noexcept
		{
			return parameters[parameterIndex];
		}

		void SetVariablesUpdated() noexcept
		{
			isActualOnVariables = false;
		}

		void SetVariable(size_t variableIndex, size_t gridIndex, FieldType value) noexcept
		{
			variables[variableIndex][gridIndex] = value;
			SetVariablesUpdated();
		}

		void SetVariable(size_t variableIndex, const Array<FieldType>& value) noexcept
		{
			Copy(value, variables[variableIndex]);
			SetVariablesUpdated();
		}

		template<typename ArrayType> requires (std::same_as<std::remove_cvref_t<ArrayType>, Array<FieldType>>)
			void SetVariables(ArrayType&& value) noexcept
		{
			variables = std::forward<ArrayType>(value);
			SetVariablesUpdated();
		}

		[[nodiscard]] const auto& GetVariable(size_t variableIndex) const noexcept
		{
			return variables[variableIndex];
		}

		[[nodiscard]] TwoLevelArray<FieldType>& GetVariables() noexcept
		{
			return variables;
		}

		[[nodiscard]] const TwoLevelArray<FieldType>& GetEquations() noexcept
		{
			Actualize();
			return equations;
		}

		void AddLocalOutputExpression(std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> expression, const std::string& name = std::string{}) noexcept
		{
			localOutputExpressions.push_back(std::make_pair(expression, name));
		}

		void AddGlobalOutputExpression(std::function<FieldType(const CurrentGlobalValues&)> expression, const std::string& name = std::string{}) noexcept
		{
			globalOutputExpressions.push_back(std::make_pair(expression, name));
		}

		void AddReductionOutputExpression
			( const std::function<FieldType(FieldType)>& externalFunction
			, const std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>& internalFunction
			, const std::string& name = std::string{}
			) noexcept
		{
			CurrentReduction reduction;
			reduction.ExternalFunction = externalFunction;
			reduction.InternalFunction = internalFunction;
			reductionOutputExpressions.push_back(std::make_pair(reduction, name));
		}

		void AddIntegralOutputExpression
			( const std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)>& integrand
			, const std::string& name = std::string{}
			) noexcept
		{
			CurrentReduction reduction;
			reduction.InternalFunction =
				[=](const auto& locals, const auto& globals) {return locals.IntegrationWeight * integrand(locals, globals); };
			reductionOutputExpressions.push_back(std::make_pair(reduction, name));
		}

		void AddPointOutputExpression(const std::array<CoordinateType, Dimension>& point, std::function<FieldType(const CurrentLocalValues&, const CurrentGlobalValues&)> expression, const std::string& name = std::string{}) noexcept
		{
			pointOutputExpressions.emplace_back(name, expression, point, discretizer->GetInterpolationWeightsVector(*grid, point));
		}

		[[nodiscard]] size_t LocalOutputExpressionCount() const noexcept
		{
			return localOutputExpressions.size();
		}

		[[nodiscard]] size_t GlobalOutputExpressionCount() const noexcept
		{
			return globalOutputExpressions.size();
		}

		[[nodiscard]] size_t ReductionOutputExpressionCount() const noexcept
		{
			return reductionOutputExpressions.size();
		}

		[[nodiscard]] size_t PointOutputExpressionCount() const noexcept
		{
			return pointOutputExpressions.size();
		}

		void CalculateLocalOutput(size_t index, Array<FieldType>& output) noexcept
		{
			Actualize();
			const auto globals = ConstructGlobalValues();
#pragma omp parallel 
			{
				auto locals = ConstructLocalValues();
#pragma omp for
				for (size_t i = 0; i < grid->GetSize(); i++)
				{
					FillAllLocals(i, locals);
					output[i] = localOutputExpressions[index].first(locals, globals);
				}
			}
		}

		[[nodiscard]] FieldType CalculateGlobalOutput(size_t index) noexcept
		{
			Actualize();
			const auto globals = ConstructGlobalValues();
			return globalOutputExpressions[index](globals);
		}

		[[nodiscard]] FieldType CalculateReductionOutput(size_t index) noexcept
		{
			Actualize();
			auto locals = ConstructLocalValues();
			const auto globals = ConstructGlobalValues();
			FieldType result = 0;
			for (size_t i = 0; i < grid->GetSize(); i++)
			{
				FillAllLocals(i, locals);
				result += reductionOutputExpressions[index].first.InternalFunction(locals, globals);
			}
			return reductionOutputExpressions[index].first.ExternalFunction(result);
		}

		[[nodiscard]] FieldType CalculatePointOutput(size_t index) noexcept
		{
			Actualize();
			auto locals = ConstructLocalValues();
			const auto globals = ConstructGlobalValues();
			const auto coefficients = grid->GetInterpolationCoefficients(pointOutputExpressions[index].point);
			FieldType result = 0;
			for (size_t i = 0; i < coefficients.size(); ++i)
			{
				FillAllLocals(coefficients[i].first, locals);
				result += pointOutputExpressions[index].expression(locals, globals) * coefficients[i].second;
			}
			return result;
		}

		void CalculateLocalOutput(TwoLevelArray<FieldType>& output) noexcept
		{
			Actualize();
			const auto globals = ConstructGlobalValues();
#pragma omp parallel
			{
				auto locals = ConstructLocalValues();
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); i++)
				{
					FillAllLocals(i, locals);
					for (size_t j = 0; j < LocalOutputExpressionCount(); j++)
					{
						output[j][i] = localOutputExpressions[j].first(locals, globals);
					}
				}
			}
		}

		void CalculateGlobalOutput(Array<FieldType>& output) noexcept
		{
			Actualize();
			const auto globals = ConstructGlobalValues();
			for (size_t i = 0; i < GlobalOutputExpressionCount(); i++)
			{
				output[i] = globalOutputExpressions[i].first(globals);
			}
		}

		void CalculateReductionOutput(Array<FieldType>& output) noexcept
		{
			Actualize();
			auto locals = ConstructLocalValues();
			const auto globals = ConstructGlobalValues();
			std::fill_n(output.begin(), ReductionOutputExpressionCount(), 0);
			for (size_t i = 0; i < grid->GetSize(); i++)
			{
				FillAllLocals(i, locals);
				for (size_t j = 0; j < ReductionOutputExpressionCount(); j++)
				{
					output[j] += reductionOutputExpressions[j].first.InternalFunction(locals, globals);
				}
			}
			for (size_t j = 0; j < ReductionOutputExpressionCount(); j++)
			{
				output[j] = reductionOutputExpressions[j].first.ExternalFunction(output[j]);
			}
		}

		void CalculatePointOutput(Array<FieldType>& output) noexcept
		{
			Actualize();
			auto locals = ConstructLocalValues();
			const auto globals = ConstructGlobalValues();
			for (size_t expressionIndex = 0; expressionIndex < PointOutputExpressionCount(); ++expressionIndex)
			{
				const auto& weights = pointOutputExpressions[expressionIndex].weights;
				for (size_t i = 0; i < weights.NonZeroCount(); ++i)
				{
					FillAllLocals(weights.GetIndex(i), locals);
					output[expressionIndex] += weights.GetValue(i) * pointOutputExpressions[expressionIndex].expression(locals, globals);
				}
			}
		}

		[[nodiscard]] std::string MakeStateName() const noexcept
		{
			std::string result = descriptor.GetProblemName();
			for (size_t parameterIndex = 0; parameterIndex < ParameterCount(); ++parameterIndex)
			{
				result.append(Format("_{}={}", descriptor.GetParameterName(parameterIndex), parameters[parameterIndex]));
			}
			return result;
		}

		static constexpr Serializer::ProblemType SerializerProblemType = Serializer::ProblemType::StationaryProblem;

		struct SerializedData
		{
			uint32_t type = static_cast<uint32_t>(SerializerProblemType);
			struct Header
			{
				uint64_t dataType = static_cast<uint32_t>(Serializer::ToDataType<FieldType>());
				uint64_t fieldCount;
				uint64_t variableCount;
				uint64_t parameterCount;
				uint64_t localOutputCount;
				uint64_t globalOutputCount;
				uint64_t reductionOutputCount;
				uint64_t pointOutputCount;

				Header() noexcept = default;

				Header(const Self& problem) noexcept
					: fieldCount(problem.GetDescriptor().ContinuousEquationCount())
					, variableCount(problem.GetDescriptor().DiscreteEquationCount())
					, parameterCount(problem.GetDescriptor().ParameterCount())
					, localOutputCount(problem.LocalOutputExpressionCount())
					, globalOutputCount(problem.GlobalOutputExpressionCount())
					, reductionOutputCount(problem.ReductionOutputExpressionCount())
					, pointOutputCount(problem.PointOutputExpressionCount())
				{}
			} header;
			Array<uint8_t> data;
		};

	protected:
		[[nodiscard]] size_t GetSerializedDataStringLength() const noexcept
		{
			size_t stringLength = 0;
			for (size_t index = 0; index < descriptor.EquationCount(); ++index)
			{
				stringLength += descriptor.GetVariableName(index).size() + 1;
			}
			for (size_t index = 0; index < descriptor.ParameterCount(); ++index)
			{
				stringLength += descriptor.GetParameterName(index).size() + 1;
			}
			for (size_t index = 0; index < LocalOutputExpressionCount(); ++index)
			{
				stringLength += localOutputExpressions[index].second.size() + 1;
			}
			for (size_t index = 0; index < GlobalOutputExpressionCount(); ++index)
			{
				stringLength += globalOutputExpressions[index].second.size() + 1;
			}
			for (size_t index = 0; index < ReductionOutputExpressionCount(); ++index)
			{
				stringLength += reductionOutputExpressions[index].second.size() + 1;
			}
			for (size_t index = 0; index < PointOutputExpressionCount(); ++index)
			{
				stringLength += pointOutputExpressions[index].name.size() + 1;
			}
			return stringLength;
		}

		[[nodiscard]] size_t GetSerializedDataFieldCount() const noexcept
		{
			return descriptor.ContinuousEquationCount() + LocalOutputExpressionCount();
		}

		[[nodiscard]] size_t GetSerializedDataVariableCount() const noexcept
		{
			return descriptor.DiscreteEquationCount() + GlobalOutputExpressionCount()
				+ ReductionOutputExpressionCount() + PointOutputExpressionCount();
		}

		uint8_t* WriteSerializedDataFields(uint8_t* data) const noexcept
		{
			auto diff = sizeof(FieldType) * variables.Flatten().size();
			std::memcpy(data, variables.Flatten().data(), diff);
			return data + diff;
		}

		uint8_t* WriteSerializedDataLocalOutput(uint8_t* data) noexcept
		{
			TwoLevelArray<FieldType> localOutput({ {LocalOutputExpressionCount(), grid->GetSize()} });
			CalculateLocalOutput(localOutput);
			auto diff = sizeof(FieldType) * localOutput.Flatten().size();
			std::memcpy(data, localOutput.Flatten().data(), diff);
			return data + diff;
		}

		uint8_t* WriteSerializedDataParameters(uint8_t* data) const noexcept
		{
			auto diff = sizeof(FieldType) * descriptor.ParameterCount();
			std::memcpy(data, parameters.data(), diff);
			return data + diff;
		}

		uint8_t* WriteSerializedDataGlobalOutput(uint8_t* data) noexcept
		{
			Array<FieldType> globalOutput(GlobalOutputExpressionCount());
			CalculateGlobalOutput(globalOutput);
			auto diff = sizeof(FieldType) * GlobalOutputExpressionCount();
			std::memcpy(data, globalOutput.data(), diff);
			data += diff;

			Array<FieldType> reductionOutput(ReductionOutputExpressionCount());
			CalculateReductionOutput(reductionOutput);
			diff = sizeof(FieldType) * ReductionOutputExpressionCount();
			std::memcpy(data, reductionOutput.data(), diff);
			data += diff;

			Array<FieldType> pointOutput(PointOutputExpressionCount());
			CalculatePointOutput(pointOutput);
			diff = sizeof(FieldType) * PointOutputExpressionCount();
			std::memcpy(data, pointOutput.data(), diff);
			return data + diff;
		}

		uint8_t* WriteSerializedDataStrings(uint8_t* data) const noexcept
		{
			size_t diff = 0;
			for (size_t index = 0; index < descriptor.EquationCount(); ++index)
			{
				diff = descriptor.GetVariableName(index).size() + 1;
				std::memcpy(data, descriptor.GetVariableName(index).c_str(), diff);
				data += diff;
			}
			for (size_t index = 0; index < descriptor.ParameterCount(); ++index)
			{
				diff = descriptor.GetParameterName(index).size() + 1;
				std::memcpy(data, descriptor.GetParameterName(index).c_str(), diff);
				data += diff;
			}
			for (size_t index = 0; index < LocalOutputExpressionCount(); ++index)
			{
				diff = localOutputExpressions[index].second.size() + 1;
				std::memcpy(data, localOutputExpressions[index].second.c_str(), diff);
				data += diff;
			}
			for (size_t index = 0; index < GlobalOutputExpressionCount(); ++index)
			{
				diff = globalOutputExpressions[index].second.size() + 1;
				std::memcpy(data, globalOutputExpressions[index].second.c_str(), diff);
				data += diff;
			}
			for (size_t index = 0; index < ReductionOutputExpressionCount(); ++index)
			{
				diff = reductionOutputExpressions[index].second.size() + 1;
				std::memcpy(data, reductionOutputExpressions[index].second.c_str(), diff);
				data += diff;
			}
			for (size_t index = 0; index < PointOutputExpressionCount(); ++index)
			{
				diff = pointOutputExpressions[index].name.size() + 1;
				std::memcpy(data, pointOutputExpressions[index].name.c_str(), diff);
				data += diff;
			}
			return data + diff;
		}

		void LoadParameters(FieldType* sourceData) noexcept
		{
			Copy(std::span(sourceData, descriptor.ParameterCount()), parameters);
		}

		void LoadVariables(FieldType* sourceData, const std::shared_ptr<Grid<Dimension, CoordinateType>>& loadedGrid) noexcept
		{
			if (*loadedGrid != *grid)
			{
				Vector<FieldType> tmpFields(loadedGrid->GetSize() * descriptor.ContinuousEquationCount() + descriptor.DiscreteEquationCount());
				Copy(std::span(sourceData, tmpFields.size()), tmpFields);
				ProblemUtils::RebaseFields(tmpFields, *loadedGrid, *this);
			}
			else
			{
				Copy(std::span(sourceData, DOFCount()), variables.Flatten());
			}
		}

	public:
		[[nodiscard]] SerializedData Save() noexcept
		{
			SerializedData result;
			result.header = SerializedData::Header(*this);

			const auto dataSize = sizeof(FieldType) * 
				(descriptor.ParameterCount() + grid->GetSize() * GetSerializedDataFieldCount() + GetSerializedDataVariableCount())
				+ GetSerializedDataStringLength();

			result.data = Array<uint8_t>(dataSize);
			auto data = result.data.data();
			data = WriteSerializedDataParameters(data);
			data = WriteSerializedDataFields(data);
			data = WriteSerializedDataLocalOutput(data);
			data = WriteSerializedDataGlobalOutput(data);
			WriteSerializedDataStrings(data);

			return result;
		}

		void Load(const Array<uint8_t>& headerData, const Array<uint8_t>& data, const std::shared_ptr<Grid<Dimension, CoordinateType>>& loadedGrid, Serializer::DataToLoad dataToLoad) noexcept
		{
			AssertE(headerData.size() == sizeof(SerializedData::Header), MessageTag::Serialization, "Invalid deserialized problem header size!");
			const auto header = reinterpret_cast<typename SerializedData::Header*>(headerData.data());
			AssertE(header->dataType == static_cast<uint32_t>(Serializer::ToDataType<FieldType>()), MessageTag::Serialization, "Invalid deserialized problem data type!");
			AssertE(header->fieldCount == descriptor.ContinuousEquationCount(), MessageTag::Serialization, "Invalid deserialized problem field count!");
			AssertE(header->variableCount == descriptor.DiscreteEquationCount(), MessageTag::Serialization, "Invalid deserialized problem variable count!");
			AssertE(header->parameterCount == descriptor.ParameterCount(), MessageTag::Serialization, "Invalid deserialized problem parameter count!");
			auto sourceData = reinterpret_cast<FieldType*>(data.data());
			if (dataToLoad & Serializer::DataToLoad::Parameters)
			{
				LoadParameters(sourceData);
			}
			sourceData += descriptor.ParameterCount();
			if (dataToLoad & Serializer::DataToLoad::Variables)
			{
				LoadVariables(sourceData, loadedGrid);
			}
		}

		~BaseProblem() = default;
	};
}