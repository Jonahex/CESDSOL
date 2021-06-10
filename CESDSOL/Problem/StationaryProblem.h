#pragma once

#include "Problem/BaseProblem.h"
#include "Problem/StationaryProblemDescriptor.h"

namespace CESDSOL
{
	template<size_t DimensionArg, template<typename> typename MatrixTypeArg = CSRMatrix, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class StationaryProblem
		: public BaseProblem<StationaryProblemDescriptor<DimensionArg, CoordinateTypeArg, FieldTypeArg>, MatrixTypeArg>
	{
	private:
		using BaseType = BaseProblem<StationaryProblemDescriptor<DimensionArg, CoordinateTypeArg, FieldTypeArg>, MatrixTypeArg>;
		
	public:
		using typename BaseType::DescriptorType;
		using typename BaseType::DiscretizationType;
		using typename BaseType::GridType;
		using typename BaseType::CoordinateType;
		using typename BaseType::FieldType;
		using typename BaseType::CurrentGlobalValuesForVIEs;

		using CurrentLocalValuesForJacobian = typename DescriptorType::CurrentLocalValuesForJacobian;
		using CurrentGlobalValuesForJacobian = typename DescriptorType::CurrentGlobalValuesForJacobian;

		using JacobianMatrixType = MatrixTypeArg<FieldType>;

		using BaseType::DOFCount;

	private:
		using BaseType::descriptor;
		using BaseType::grid;

		using BaseType::differentiationWeights;
		using BaseType::fieldDerivativeOperatorMap;

		using BaseType::parameters;
		using BaseType::variables;
		using BaseType::equations;
		using BaseType::globalPIEs;
		using BaseType::globalVIEs;
		using BaseType::globalVDEs;
		using BaseType::reductions;

		using BaseType::Actualize;
		using BaseType::ConstructDerivativesLevelStructure;
		using BaseType::GetTrueRegionIndex;
		using BaseType::FillAllLocals;
		
		friend DescriptorType;

		bool isJacobianReady = false;
		ThreeLevelArray<Array<FieldType>> jacobian;
		ThreeLevelArray<Array<FieldType>> lvdeJacobians;
		TwoLevelArray<FieldType> gvdeJacobians;
		ThreeLevelArray<Array<FieldType>> reductionJacobians;

		struct JacobianElement
		{
			size_t Index;
			size_t OperatorIndex;
			CoordinateType Multiplier;
		};

		ThreeLevelArray<Array<JacobianElement>> jacobianStructure;
		JacobianMatrixType jacobianMatrix;
		bool isActualOnJacobian = false;

		[[nodiscard]] ThreeLevelArray<Array<FieldType>> ConstructJacobianCommon(size_t size) const noexcept
		{
			Array<std::pair<size_t, LevelStructure<2>>> levelStructure{ size };
			for (size_t i = 0; i < size; i++)
			{
				Array<std::pair<size_t, size_t>> internalStructure{ descriptor.EquationCount() };
				for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
				{
					internalStructure[j] = { 1, descriptor.DerivativeOperatorCount(j) + 1 };
				}
				for (size_t j = 0; j < descriptor.DiscreteEquationCount(); j++)
				{
					internalStructure[descriptor.ContinuousEquationCount() + j] = { 1, 1 };
				}
				levelStructure[i] = { 1, internalStructure };
			}
			auto result = ThreeLevelArray<Array<FieldType>>(levelStructure);
			for (size_t i = 0; i < size; i++)
			{
				for (size_t j = 0; j < descriptor.EquationCount(); j++)
				{
					for (size_t k = 0; k <= (j < descriptor.ContinuousEquationCount() ? descriptor.DerivativeOperatorCount(j) : 0); k++)
					{
						result[i][j][k] = Array<FieldType>(grid->GetSize());
					}
				}
			}
			return result;
		}

		[[nodiscard]] ThreeLevelArray<Array<FieldType>> ConstructJacobian() const noexcept
		{
			return ConstructJacobianCommon(descriptor.EquationCount());
		}

		void FillJacobian() noexcept
		{
			for (size_t i = 0; i < descriptor.EquationCount(); i++)
			{
				for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
				{
					for (size_t k = 0; k <= descriptor.DerivativeOperatorCount(j); k++)
					{
						for (size_t l = 0; l < grid->GetRegionCount(); l++)
						{
							if (descriptor.HasJacobianComponent(i, j, k, l))
							{
								jacobian[i][j][k] = Array<FieldType>(grid->GetSize());
								break;
							}
						}
					}
				}
				for (size_t j = 0; j < descriptor.DiscreteEquationCount(); j++)
				{
					if (descriptor.HasJacobianComponent(i, j, 0, 0))
					{
						jacobian[i][j][0] = Array<FieldType>(1);
					}
				}
			}
		}

		[[nodiscard]] ThreeLevelArray<Array<FieldType>> ConstructReductionJacobian() noexcept
		{
			return ConstructJacobianCommon(descriptor.ReductionCount());
		}

		void FillReductionJacobian() noexcept
		{
			for (size_t i = 0; i < descriptor.ReductionCount(); i++)
			{
				for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
				{
					for (size_t k = 0; k <= descriptor.DerivativeOperatorCount(j); k++)
					{
						if (descriptor.HasReductionJacobianComponent(i, j, k))
						{
							reductionJacobians[i][j][k] = Array<FieldType>(grid->GetSize());
						}
					}
				}
				for (size_t j = 0; j < descriptor.DiscreteEquationCount(); j++)
				{
					if (descriptor.HasReductionJacobianComponent(i, j, 0))
					{
						reductionJacobians[i][j][0] = Array<FieldType>(1);
					}
				}
			}
		}

		[[nodiscard]] ThreeLevelArray<Array<FieldType>> ConstructLVDEJacobian() noexcept
		{
			return ConstructJacobianCommon(descriptor.LocalVDECount());
		}

		void FillLVDEJacobian() noexcept
		{
			for (size_t i = 0; i < descriptor.LocalVDECount(); i++)
			{
				for (size_t j = 0; j < descriptor.ContinuousEquationCount(); j++)
				{
					for (size_t k = 0; k <= j < descriptor.ContinuousEquationCount() ? descriptor.DerivativeOperatorCount(j) : 0; k++)
					{
						if (descriptor.HasLVDEJacobianComponent(i, j, k))
						{
							lvdeJacobians[i][j][k] = Array<FieldType>(grid->GetSize());
						}
					}
				}
			}
		}

		CurrentLocalValuesForJacobian ConstructLocalValuesForJacobian() const noexcept
		{
			return CurrentLocalValuesForJacobian
			{
				descriptor.LocalPIECount(),
				descriptor.LocalVIECount(),
				descriptor.ContinuousEquationCount(),
				ConstructDerivativesLevelStructure(),
				descriptor.LocalVDECount(),
				ApproximateStructure(lvdeJacobians),
				ApproximateStructure(reductionJacobians)
			};
		}

		CurrentGlobalValuesForJacobian ConstructGlobalValuesForJacobian() const noexcept
		{
			return CurrentGlobalValuesForJacobian
			{
				globalPIEs,
				parameters,
				globalVIEs,
				std::span(variables[descriptor.ContinuousEquationCount()].data(), variables[descriptor.ContinuousEquationCount()].data() + descriptor.DiscreteEquationCount()),
				globalVDEs,
				reductions,
				gvdeJacobians
			};
		}

		void FillLVDEJacobians(size_t pointIndex, const CurrentLocalValuesForJacobian& locals) noexcept
		{
			for (size_t j = 0; j < descriptor.LocalVDECount(); j++)
			{
				for (size_t k = 0; k < descriptor.EquationCount(); k++)
				{
					for (size_t l = 0; l <= (k < descriptor.ContinuousEquationCount() ? descriptor.DerivativeOperatorCount(k) : 0); ++l)
					{
						locals.LVDEJacobianComponentValues[j][k][l] = lvdeJacobians[j][k][l][pointIndex];
					}
				}
			}
		}

		void FillReductionJacobians(size_t pointIndex, const CurrentLocalValuesForJacobian& locals) noexcept
		{
			for (size_t j = 0; j < descriptor.ReductionCount(); j++)
			{
				for (size_t k = 0; k < descriptor.ContinuousEquationCount(); k++)
				{
					for (size_t l = 0; l <= descriptor.DerivativeOperatorCount(k); l++)
					{
						locals.ReductionJacobianComponentValues[j][k][l] = reductionJacobians[j][k][l][pointIndex];
					}
				}
			}
		}

		void FillReductionJacobiansGlobal(const CurrentLocalValuesForJacobian& locals) noexcept
		{
			for (size_t j = 0; j < descriptor.ReductionCount(); ++j)
			{
				for (size_t k = descriptor.ContinuousEquationCount(); k < descriptor.EquationCount(); ++k)
				{
					locals.ReductionJacobianComponentValues[j][k][0] = reductionJacobians[j][k][0][0];
				}
			}
		}

		void UpdateExpressionJacobians() noexcept
		{
			const auto globals = ConstructGlobalValuesForJacobian();
			for (size_t i = 0; i < descriptor.GlobalVDECount(); ++i)
			{
				for (size_t j = 0; j < descriptor.DiscreteEquationCount(); ++j)
				{
					gvdeJacobians[i][j] = descriptor.CalculateGVDEJacobianComponent(i, j, globals);
				}
			}
			for (size_t j = 0; j < descriptor.LocalVDECount(); ++j)
			{
				for (size_t k = descriptor.ContinuousEquationCount(); k < descriptor.EquationCount(); ++k)
				{
					reductionJacobians[j][k][0][0] = 0;
				}
			}
#pragma omp parallel
			{
				auto locals = ConstructLocalValuesForJacobian();
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); ++i)
				{
					FillAllLocals(i, locals);
					for (size_t j = 0; j < descriptor.LocalVDECount(); ++j)
					{
						for (size_t k = 0; k < descriptor.EquationCount(); ++k)
						{
							for (size_t l = 0; l <= (k < descriptor.ContinuousEquationCount() ? descriptor.DerivativeOperatorCount(k) : 0); ++l)
							{
								if (descriptor.HasLVDEJacobianComponent(j, k, l))
								{
									lvdeJacobians[j][k][l][i] = descriptor.CalculateLVDEJacobianComponent(j, k, l, locals, globals);
									locals.LVDEJacobianComponentValues[j][k][l] = lvdeJacobians[j][k][l][i];
								}
							}
						}
					}
					for (size_t j = 0; j < descriptor.ReductionCount(); j++)
					{
						for (size_t k = 0; k < descriptor.ContinuousEquationCount(); k++)
						{
							if (descriptor.HasReductionJacobianComponent(j, k, 0))
							{
								reductionJacobians[j][k][0][i] = descriptor.CalculateReductionJacobianComponent(j, k, 0, locals, globals);
							}
							for (size_t l = 1; l <= descriptor.DerivativeOperatorCount(k); l++)
							{
								if (descriptor.HasReductionJacobianComponent(j, k, l))
								{
									const auto& weightsMatrix = differentiationWeights[fieldDerivativeOperatorMap[k][l - 1]];
									for (size_t n = weightsMatrix.GetRowCount(i); n < weightsMatrix.GetRowCount(i + 1); ++n)
									{
										reductionJacobians[j][k][l][weightsMatrix.GetColumnIndex(n)] +=
											descriptor.CalculateReductionJacobianComponent(j, k, l, locals, globals) * weightsMatrix.GetValue(n);
									}
								}
							}
						}
						for (size_t k = descriptor.ContinuousEquationCount(); k < descriptor.EquationCount(); ++k)
						{
							if (descriptor.HasReductionJacobianComponent(j, k, 0))
							{
								reductionJacobians[j][k][0][0] += descriptor.CalculateReductionInternalJacobianComponent(j, k, 0, locals, globals);
							}
						}
					}
				}
			}
			for (size_t j = 0; j < descriptor.LocalVDECount(); ++j)
			{
				for (size_t k = descriptor.ContinuousEquationCount(); k < descriptor.EquationCount(); ++k)
				{
					reductionJacobians[j][k][0][0] *= descriptor.CalculateReductionExternalJacobianComponent(j, globals);
				}
			}
		}

		void CalculateJacobian() noexcept
		{
			const auto globals = ConstructGlobalValuesForJacobian();
#pragma omp parallel
			{
				auto locals = ConstructLocalValuesForJacobian();
				FillReductionJacobiansGlobal(locals);
#pragma omp for
				for (int64_t i = 0; i < grid->GetSize(); i++)
				{
					FillAllLocals(i, locals);
					FillLVDEJacobians(i, locals);
					FillReductionJacobians(i, locals);
					const auto regionIndex = grid->GetRegionIndex(i);
					for (size_t j = 0; j < descriptor.EquationCount(); j++)
					{
						auto trueRegionIndex = regionIndex;
						if (j >= descriptor.ContinuousEquationCount() || !descriptor.HasContinuousEquation(j, regionIndex))
						{
							trueRegionIndex = 0;
						}
						for (size_t k = 0; k < descriptor.EquationCount(); k++)
						{
							for (size_t l = 0; l <= (k < descriptor.ContinuousEquationCount() ? descriptor.DerivativeOperatorCount(k) : 0); l++)
							{
								if (descriptor.HasJacobianComponent(j, k, l, trueRegionIndex))
								{
									jacobian[j][k][l][i] = descriptor.CalculateJacobianComponent(j, k, l, trueRegionIndex, locals, globals);
								}
							}
						}
					}
				}
			}
		}

		void CalculateJacobianStructure() noexcept
		{
			if constexpr (std::is_same_v<JacobianMatrixType, CSRMatrix<CoordinateType>>)
			{
				const auto dofCount = DOFCount();
				const auto ceCount = descriptor.ContinuousEquationCount();
				const auto deCount = descriptor.DiscreteEquationCount();
				const auto eCount = descriptor.EquationCount();

				jacobianStructure = ThreeLevelArray<Array<JacobianElement>>(
					{ {ceCount, {grid->GetSize(), ceCount}} });

				std::atomic<size_t> nonzeroCount = 0;
				for (size_t i = 0; i < ceCount; ++i)
				{
#pragma omp parallel for
					for (int64_t k = 0; k < grid->GetSize(); ++k)
					{
						const auto trueRegionIndex = GetTrueRegionIndex(i, k);
						for (size_t j = 0; j < ceCount; ++j)
						{
							size_t elementCount = 0;
							if (descriptor.HasJacobianComponent(i, j, 0, trueRegionIndex))
							{
								++elementCount;
							}
							for (size_t l = 1; l <= descriptor.DerivativeOperatorCount(j); ++l)
							{
								if (descriptor.HasJacobianComponent(i, j, l, trueRegionIndex))
								{
									elementCount += differentiationWeights[fieldDerivativeOperatorMap[j][l - 1]].GetRowLength(k);
								}
							}
							auto& row = jacobianStructure[i][k][j] = Array<JacobianElement>(elementCount);
							size_t elementIndex = 0;
							if (descriptor.HasJacobianComponent(i, j, 0, trueRegionIndex))
							{
								row[elementIndex++] = { j * grid->GetSize() + k, 0, 1 };
							}
							for (size_t l = 1; l <= descriptor.DerivativeOperatorCount(j); ++l)
							{
								if (descriptor.HasJacobianComponent(i, j, l, trueRegionIndex))
								{
									const auto& weightMatrix = differentiationWeights[fieldDerivativeOperatorMap[j][l - 1]];
									for (size_t m = weightMatrix.GetRowCount(k); m < weightMatrix.GetRowCount(k + 1); ++m)
									{
										row[elementIndex++] = { j * grid->GetSize() + weightMatrix.GetColumnIndex(m), l, weightMatrix.GetValue(m) };
									}
								}
							}
							std::sort(row.begin(), row.end(), [](const auto& left, const auto& right) {return left.Index < right.Index; });

							if (elementCount > 0)
							{
								size_t setElements = 1;
								auto currentElement = row[0].Index;
								for (size_t l = 1; l < elementCount; l++)
								{
									if (row[l].Index != currentElement)
									{
										++setElements;
										currentElement = row[l].Index;
									}
								}
								nonzeroCount += setElements;
							}
						}

						for (size_t j = ceCount; j < eCount; ++j)
						{
							if (descriptor.HasJacobianComponent(i, j, 0, trueRegionIndex))
							{
								++nonzeroCount;
							}
						}
					}
				}

				for (size_t i = ceCount; i < eCount; ++i)
				{
					for (size_t j = 0; j < ceCount; ++j)
					{
						for (size_t k = 0; k < jacobian[i][j].size(); k++)
						{
							if (descriptor.HasJacobianComponent(i, j, k, 0))
							{
								nonzeroCount += grid->GetSize();
								break;
							}
						}
					}
					for (size_t j = ceCount; j < eCount; ++j)
					{
						if (descriptor.HasJacobianComponent(i, j, 0, 0))
						{
							++nonzeroCount;
						}
					}
				}

				jacobianMatrix = CSRMatrix<FieldType>(dofCount, dofCount, nonzeroCount);
				size_t setElements = 0;
				for (size_t i = 0; i < ceCount; ++i)
				{
					for (size_t j = 0; j < grid->GetSize(); ++j)
					{
						const auto rowIndex = i * grid->GetSize() + j;
						jacobianMatrix.SetRowCount(rowIndex, setElements);
						for (size_t k = 0; k < ceCount; ++k)
						{
							auto& row = jacobianStructure[i][j][k];
							if (row.size() > 0)
							{
								auto currentElement = row[0].Index;
								row[0].Index = setElements;
								jacobianMatrix.SetColumnIndex(setElements++, currentElement);
								for (size_t l = 1; l < row.size(); l++)
								{
									if (currentElement != row[l].Index)
									{
										currentElement = row[l].Index;
										row[l].Index = setElements;
										jacobianMatrix.SetColumnIndex(setElements++, currentElement);
									}
									else
									{
										row[l].Index = setElements - 1;
									}
								}
							}
						}
						for (size_t k = ceCount; k < eCount; ++k)
						{
							const auto trueRegionIndex = GetTrueRegionIndex(i, j);
							if (descriptor.HasJacobianComponent(i, k, 0, trueRegionIndex))
							{
								jacobianMatrix.SetColumnIndex(setElements++, ceCount * grid->GetSize() + k - ceCount);
							}
						}
					}
				}
				for (size_t i = ceCount; i < eCount; ++i)
				{
					jacobianMatrix.SetRowCount(ceCount * grid->GetSize() + i - ceCount, setElements);
					for (size_t j = 0; j < ceCount; ++j)
					{
						for (size_t k = 0; k < jacobian[i][j].size(); ++k)
						{
							if (descriptor.HasJacobianComponent(i, j, k, 0))
							{
								for (size_t l = 0; l < grid->GetSize(); ++l)
								{
									jacobianMatrix.SetColumnIndex(setElements++, j * grid->GetSize() + l);
								}
								break;
							}
						}
					}
					for (size_t j = ceCount; j < eCount; ++j)
					{
						if (descriptor.HasJacobianComponent(i, j, 0, 0))
						{
							jacobianMatrix.SetColumnIndex(setElements++, ceCount * grid->GetSize() + j - ceCount);
						}
					}
				}
				AssertE(setElements == nonzeroCount, MessageTag::Problem, "Error in jacobian structure calculation!");
				jacobianMatrix.SetRowCount(dofCount, setElements);
			}
		}

		void UpdateJacobian() noexcept
		{
			UpdateExpressionJacobians();
			CalculateJacobian();

			if constexpr (std::is_same_v<JacobianMatrixType, CSRMatrix<CoordinateType>>)
			{
				const auto ceCount = descriptor.ContinuousEquationCount();
				const auto eCount = descriptor.EquationCount();

				jacobianMatrix.Nullify();
				for (size_t i = 0; i < ceCount; i++)
				{
#pragma omp parallel for
					for (int64_t j = 0; j < grid->GetSize(); j++)
					{
						auto lastIndex = jacobianMatrix.GetRowCount(i * grid->GetSize() + j);
						for (size_t k = 0; k < ceCount; k++)
						{
							const auto& row = jacobianStructure[i][j][k];
							for (size_t l = 0; l < row.size(); l++)
							{
								lastIndex = row[l].Index;
								jacobianMatrix.SetValue(lastIndex, jacobianMatrix.GetValue(lastIndex) + row[l].Multiplier * jacobian[i][k][row[l].OperatorIndex][j]);
							}
						}
						for (size_t k = ceCount; k < eCount; ++k)
						{
							const auto trueRegionIndex = GetTrueRegionIndex(i, j);
							if (descriptor.HasJacobianComponent(i, k, 0, trueRegionIndex))
							{
								jacobianMatrix.SetValue(lastIndex++, jacobian[i][k][0][j]);
							}
						}
					}
				}
				for (size_t i = ceCount; i < eCount; ++i)
				{
					auto lastIndex = jacobianMatrix.GetRowCount(ceCount * grid->GetSize() + i - ceCount);
					for (size_t k = 0; k < ceCount; k++)
					{
						if (descriptor.HasJacobianComponent(i, k, 0, 0))
						{
							for (size_t l = 0; l < grid->GetSize(); ++l)
							{
								jacobianMatrix.SetValue(lastIndex++, jacobian[i][k][0][l]);
							}
						}
					}
					for (size_t k = ceCount; k < eCount; ++k)
					{
						if (descriptor.HasJacobianComponent(i, k, 0, 0))
						{
							jacobianMatrix.SetValue(lastIndex++, jacobian[i][k][0][0]);
						}
					}
				}
			}
		}

		StationaryProblem(sptr<GridType> grid, uptr<DiscretizationType> discretizer, const DescriptorType& descriptor)
		: BaseType(std::move(grid), std::move(discretizer), descriptor)
		, jacobian(ConstructJacobian())
		, lvdeJacobians(ConstructLVDEJacobian())
		, gvdeJacobians({ {descriptor.LocalVDECount(), descriptor.DiscreteEquationCount()} })
		, reductionJacobians(ConstructReductionJacobian())
		{}

	public:
		[[nodiscard]] FieldType CalculateSolutionNorm() noexcept
		{
			Actualize();
			return descriptor.CalculateMerit(variables.Flatten());
		}

		[[nodiscard]] FieldType GetMerit() noexcept
		{
			Actualize();
			return descriptor.CalculateMerit(equations.Flatten());
		}
		
		[[nodiscard]] const JacobianMatrixType& GetJacobian() noexcept
		{
			Actualize();
			if (!isActualOnJacobian)
			{
				if (jacobianStructure.size() == 0)
				{
					CalculateJacobianStructure();
				}
				UpdateJacobian();
			}
			return jacobianMatrix;
		}

#ifdef DebugMode
		void PrintJacobianStructure(std::ostream& stream) noexcept
		{
			if (jacobianStructure.size() == 0)
			{
				CalculateJacobianStructure();
			}
			for (size_t i = 0; i < jacobianStructure.size(); ++i)
			{
				for (size_t j = 0; j < jacobianStructure[i].size(); ++j)
				{
					for (size_t k = 0; k < jacobianStructure[i][j].size(); ++k)
					{
						for (size_t l = 0; l < jacobianStructure[i][j][k].size(); ++l)
						{
							const auto& element = jacobianStructure[i][j][k][l];
							stream << Format("{} {} {} {}: {} {} {}\n", i, j, k, l, element.Index, element.OperatorIndex, element.Multiplier);
						}
					}
				}
			}
		}
#endif
	};

	template<size_t Dimension, typename CoordinateType, typename FieldType>
	template<template<typename> typename MatrixType>
	[[nodiscard]] auto	
		StationaryProblemDescriptor<Dimension, CoordinateType, FieldType>::MakeProblem
		( sptr<Grid<Dimension, CoordinateType>> grid
		, uptr<Discretization<Dimension, MatrixType, CoordinateType>> discretizer
		) const noexcept
	{
		AssertE(grid->GetDescriptor() == this->GetGridDescriptor() || this->Validate(),
			MessageTag::Problem, "Invalid problem descriptor.");
		return sptr<StationaryProblem<Dimension, MatrixType, CoordinateType, FieldType>>(
			new StationaryProblem<Dimension, MatrixType, CoordinateType, FieldType>(std::move(grid), std::move(discretizer), *this));
	}
}