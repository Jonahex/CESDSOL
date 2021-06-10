#pragma once

#include "Discretization/Discretization.h"
#include "Grid/DirectProductGrid.h"
#include "Math/Concepts.h"
#include "Math/LinearAlgebra.h"
#include "Math/MultiLevelArray.h"
#include "Math/CSRMatrixOperations.h"

#include <span>

namespace CESDSOL
{
	class StructuredFiniteDifferenceDiscretizationCalculator
	{
	private:
		template<typename ScalarType, Concepts::Vector InputCollection, Concepts::Vector OutputCollection>
		constexpr void GenerateFornbergWeights(const InputCollection& grid,
			OutputCollection& weights, size_t derivativeOrder, ScalarType center = 0) const noexcept
		{
			const auto stencilSize = grid.size();
			auto tmp = TwoLevelArray<ScalarType>({ derivativeOrder + 1, stencilSize });
			tmp[0][0] = 1;
			ScalarType previousDifferenceProduct = 1;
			auto shift = grid[0] - center;
			for (size_t i = 1; i < stencilSize; ++i)
			{
				const auto mn = std::min(i, derivativeOrder);
				ScalarType differenceProduct = 1;
				const auto previousShift = shift;
				shift = grid[i] - center;
				for (size_t j = 0; j < i; ++j)
				{
					const auto difference = grid[i] - grid[j];
					differenceProduct *= difference;
					if (j == i - 1)
					{
						const auto multiplier = previousDifferenceProduct / differenceProduct;
						for (i64 k = mn; k >= 1; --k)
						{
							tmp[k][i] = multiplier * (k * tmp[k - 1][i - 1] - previousShift * tmp[k][i - 1]);
						};
						tmp[0][i] = -multiplier * previousShift * tmp[0][i - 1];
					};

					const auto invDifference = 1 / difference;
					for (i64 k = mn; k >= 1; --k)
					{
						tmp[k][j] = (shift * tmp[k][j] - k * tmp[k - 1][j]) * invDifference;
					};
					tmp[0][j] *= shift * invDifference;
				};
				previousDifferenceProduct = differenceProduct;
			};
			for (size_t i = 0; i < stencilSize; i++)
			{
				weights[i] = tmp[derivativeOrder][i];
			};
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] constexpr CSRMatrix<ScalarType> GetDifferentiationMatrix(const DirectProductGrid<Dimension, ScalarType>& grid, size_t dimensionIndex,
			size_t derivativeOrder, size_t stencilSize) const noexcept
		{
			const auto dimensionSize = grid.GetDimensionSize(dimensionIndex);
			const auto isPeriodic = grid.IsPeriodicDimension(dimensionIndex);
			const auto& gridData = grid.GetGrid(dimensionIndex);
			const auto halfStencil = stencilSize / 2;
			const auto period = grid.GetPeriod(dimensionIndex);

			auto weights = TwoLevelArray<ScalarType>({ dimensionSize, stencilSize });

			if (!isPeriodic)
			{
				for (size_t i = 0; i < halfStencil; ++i)
				{
					GenerateFornbergWeights(std::span(gridData.data(), stencilSize),
						weights[i], derivativeOrder, gridData[i]);
				}
				for (size_t i = dimensionSize - halfStencil; i < dimensionSize; ++i)
				{
					GenerateFornbergWeights(std::span(gridData.end() - stencilSize, gridData.end()), 
						weights[i], derivativeOrder, gridData[i]);
				}
			}
			else
			{
				auto tmp = Array<ScalarType>(stencilSize);
				for (size_t i = 0; i < halfStencil; i++)
				{
					for (size_t j = 0; j < halfStencil - i; j++)
					{
						tmp[j] = gridData[dimensionSize - halfStencil + i + j] - *period;
					}
					for (size_t j = halfStencil - i; j < stencilSize; j++)
					{
						tmp[j] = gridData[j - halfStencil + i];
					}
					GenerateFornbergWeights(tmp, weights[i], derivativeOrder, gridData[i]);
				}
				for (size_t i = dimensionSize - halfStencil; i < dimensionSize; i++)
				{
					for (size_t j = 0; j < dimensionSize - i + halfStencil; j++)
					{
						tmp[j] = gridData[i - halfStencil + j];
					}
					for (size_t j = dimensionSize - i + halfStencil; j < stencilSize; j++)
					{
						tmp[j] = *period + gridData[j - (dimensionSize - i + halfStencil)];
					}
					GenerateFornbergWeights(tmp, weights[i], derivativeOrder, gridData[i]);
				}
			}
			for (size_t i = halfStencil; i < dimensionSize - halfStencil; i++)
			{
				GenerateFornbergWeights(std::span(gridData.begin() + i - halfStencil, stencilSize), 
					weights[i], derivativeOrder, gridData[i]);
			};

			const size_t nonzeroCount = std::accumulate(weights.begin(), weights.end(), 0, 
				[](auto sum, const auto& element) {return sum + CountNonZero(element); }) * 
				grid.GetSize() / dimensionSize;

			auto result = CSRMatrix<ScalarType>(grid.GetSize(), grid.GetSize(), nonzeroCount);

			std::array<size_t, Dimension> tmp;
			size_t setCount = 0;
			auto setElement = [&](auto tmpWeight, auto index)
			{
				if (std::abs(tmpWeight) > std::numeric_limits<ScalarType>::epsilon())
				{
					result.SetValue(setCount, tmpWeight);
					tmp[dimensionIndex] = index;
					result.SetColumnIndex(setCount, grid.GetSingleIndexByMultiIndex(tmp));
					++setCount;
				};
			};
			for (size_t i = 0; i < grid.GetSize(); i++)
			{
				result.SetRowCount(i, setCount);
				grid.GetMultiIndexBySingleIndex(i, tmp);
				const size_t directionCoordinate = tmp[dimensionIndex];
				if (!isPeriodic || (directionCoordinate < dimensionSize - halfStencil && directionCoordinate >= halfStencil))
				{
					auto left = (directionCoordinate < halfStencil ? 0 : (directionCoordinate >= dimensionSize - halfStencil ?
						dimensionSize - stencilSize : directionCoordinate - halfStencil));
					for (size_t j = 0; j < stencilSize; ++j)
					{
						setElement(weights[directionCoordinate][j], left + j);
					};
				}
				else if (directionCoordinate >= dimensionSize - halfStencil)
				{
					for (size_t j = dimensionSize - directionCoordinate + halfStencil; j < stencilSize; ++j)
					{
						setElement(weights[directionCoordinate][j], j - (dimensionSize - directionCoordinate + halfStencil));
					}
					for (size_t j = 0; j < dimensionSize - directionCoordinate + halfStencil; ++j)
					{
						setElement(weights[directionCoordinate][j], directionCoordinate - halfStencil + j);
					};
				}
				else
				{
					for (size_t j = halfStencil - directionCoordinate; j < stencilSize; ++j)
					{
						setElement(weights[directionCoordinate][j], j - halfStencil + directionCoordinate);
					};
					for (size_t j = 0; j < halfStencil - directionCoordinate; ++j)
					{
						setElement(weights[directionCoordinate][j], dimensionSize - halfStencil + directionCoordinate + j);
					}
				}
			};

			return result;
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] constexpr SparseVector<ScalarType> GetInterpolationWeightsVector(const DirectProductGrid<Dimension, ScalarType>& grid, size_t dimensionIndex,
			ScalarType point, size_t stencilSize) const noexcept
		{
			const auto dimensionSize = grid.GetDimensionSize(dimensionIndex);
			const auto& dimensionGrid = grid.GetGrid(dimensionIndex);
			const auto centerIndex = std::min(LowerBoundIndexBinary(dimensionGrid, point), static_cast<int64_t>(dimensionSize - 1));
			auto leftIndex = std::max(static_cast<int64_t>(centerIndex - stencilSize / 2), 0ll);
			if (leftIndex + stencilSize >= dimensionSize)
			{
				leftIndex = std::max(static_cast<int64_t>(dimensionSize - stencilSize), 0ll);
			}
			const auto rightIndex = std::min(leftIndex + stencilSize - 1, dimensionSize - 1);
			const auto actualStencilSize = rightIndex - leftIndex + 1;
			auto result = SparseVector<ScalarType>(dimensionSize, actualStencilSize);
			for (size_t i = leftIndex; i <= rightIndex; i++)
			{
				result.SetIndex(i - leftIndex, i);
			}
			GenerateFornbergWeights(std::span(dimensionGrid.begin() + leftIndex, actualStencilSize), result.GetValues(), 0, point);
			return result;
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] constexpr Vector<ScalarType> GetIntegrationWeightsVector(const DirectProductGrid<Dimension, ScalarType>& grid, size_t dimensionIndex) const noexcept
		{
			const auto dimensionSize = grid.GetDimensionSize(dimensionIndex);
			const auto& dimensionGrid = grid.GetGrid(dimensionIndex);
			auto result = Vector<ScalarType>(dimensionSize);
			for (size_t i = 1; i < dimensionSize - 1; i++)
			{
				result[i] = 0.5 * (dimensionGrid[i + 1] - dimensionGrid[i - 1]);
			}
			if (grid.IsPeriodicDimension(dimensionIndex))
			{
				const auto period = *grid.GetPeriod(dimensionIndex);
				result[0] = 0.5 * (period + dimensionGrid[1] - dimensionGrid.Back());
				result[dimensionSize - 1] = 0.5 * (dimensionGrid[0] + period - dimensionGrid[dimensionSize - 2]);
			}
			else
			{
				result[0] = 0.5 * (dimensionGrid[1] - dimensionGrid[0]);
				result[dimensionSize - 1] = 0.5 * (dimensionGrid.Back() - dimensionGrid[dimensionSize - 2]);
			}
			return result;
		}

	public:
		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] CSRMatrix<ScalarType> GetDifferentiationMatrix(const DirectProductGrid<Dimension, ScalarType>& grid,
			const std::array<size_t, Dimension>& derivativeOrders, const std::array<size_t, Dimension>& stencilSizes) const noexcept
		{
			bool started = false;
			CSRMatrix<ScalarType> tmp;
			for (size_t i = 0; i < Dimension; i++)
			{
				if (derivativeOrders[i] > 0)
				{
					if (started)
					{
						tmp = Multiply(tmp, GetDifferentiationMatrix(grid, i, derivativeOrders[i], stencilSizes[i]));
					}
					else
					{
						tmp = GetDifferentiationMatrix(grid, i, derivativeOrders[i], stencilSizes[i]);	
						started = true;
					}
				}
			}
			if (started)
			{
				return tmp;
			}
			else
			{
				Logger::Log(MessageType::Warning, MessagePriority::Medium, MessageTag::Discretization,
					"Trivial differentiation matrix is constructed");
				return MakeIdentityMatrix<CSRMatrix<ScalarType>>(grid.GetSize());
			}
		}

		template<typename ScalarType>
		[[nodiscard]] CSRMatrix<ScalarType> GetDifferentiationMatrix(const DirectProductGrid<1, ScalarType>& grid,
			size_t derivativeOrder, size_t stencilSize) const noexcept
		{
			return GetDifferentiationMatrix(grid, { {derivativeOrder} }, { {stencilSize} });
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] CSRMatrix<ScalarType> GetDifferentiationMatrix(const DirectProductGrid<Dimension, ScalarType>& grid,
			const std::array<size_t, Dimension>& derivativeOrders, size_t stencilSize) const noexcept
		{
			std::array<size_t, Dimension> stencilSizes{};
			std::fill(stencilSizes.begin(), stencilSizes.end(), stencilSize);
			return GetDifferentiationMatrix(grid, derivativeOrders, stencilSizes);
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] SparseVector<ScalarType> GetInterpolationWeightsVector(const DirectProductGrid<Dimension, ScalarType>& grid,
			const std::array<ScalarType, Dimension>& point, const std::array<size_t, Dimension>& stencilSizes) const noexcept
		{
			SparseVector<ScalarType> result = GetInterpolationWeightsVector(grid, 0, point[0], stencilSizes[0]);
			for (size_t i = 1; i < Dimension; i++)
			{
				result = DirectProductAsVector(result, GetInterpolationWeightsVector(grid, i, point[i], stencilSizes[i]));
			}
			return result;
		}

		template<typename ScalarType>
		[[nodiscard]] SparseVector<ScalarType> GetInterpolationWeightsVector(const DirectProductGrid<1, ScalarType>& grid,
			const ScalarType& point, size_t stencilSize) const noexcept
		{
			return GetInterpolationWeightsVector(grid, { {point} }, { {stencilSize} });
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] SparseVector<ScalarType> GetInterpolationWeightsVector(const DirectProductGrid<Dimension, ScalarType>& grid,
			const std::array<ScalarType, Dimension>& point, size_t stencilSize) const noexcept
		{
			std::array<size_t, Dimension> stencilSizes{};
			std::fill(stencilSizes.begin(), stencilSizes.end(), stencilSize);
			return GetInterpolationWeightsVector(grid, point, stencilSizes);
		}

		template<size_t Dimension, typename ScalarType>
		[[nodiscard]] Vector<ScalarType> GetIntegrationWeightsVector(const DirectProductGrid<Dimension, ScalarType>& grid) const noexcept
		{
			auto result = GetIntegrationWeightsVector(grid, 0);
			for (size_t i = 1; i < Dimension; i++)
			{
				result = DirectProductAsVector(result, GetIntegrationWeightsVector(grid, i));
			}
			return result;
		}
	};

	template<size_t DimensionArg, typename CoordinateTypeArg = double>
	class StructuredFiniteDifferenceDiscretization : public Discretization<DimensionArg, CSRMatrix, CoordinateTypeArg>
	{
	public:
		using BaseType = Discretization<DimensionArg, CSRMatrix, CoordinateTypeArg>;
		using typename BaseType::CoordinateType;
		using typename BaseType::MatrixType;

		static constexpr size_t Dimension = BaseType::Dimension;

	public:
		StructuredFiniteDifferenceDiscretization(size_t stencilSize) noexcept
			: stencilSizes(FillStencils(stencilSize))
		{}

		StructuredFiniteDifferenceDiscretization(const std::array<size_t, Dimension>& stencilSizes) noexcept 
			: stencilSizes(stencilSizes)
		{}

		[[nodiscard]] MatrixType GetDifferentiationMatrix(const Grid<Dimension, CoordinateType>& grid,
			const std::array<size_t, Dimension>& derivativeOrders) const noexcept override
		{
			const auto& dpGrid = dynamic_cast<const DirectProductGrid<Dimension, CoordinateType>&>(grid);
			return StructuredFiniteDifferenceDiscretizationCalculator().GetDifferentiationMatrix(dpGrid, derivativeOrders, stencilSizes);
		}

		[[nodiscard]] SparseVector<CoordinateType> GetInterpolationWeightsVector(const Grid<Dimension, CoordinateType>& grid,
			const std::array<CoordinateType, Dimension>& point) const noexcept override
		{
			const auto& dpGrid = dynamic_cast<const DirectProductGrid<Dimension, CoordinateType>&>(grid);
			return StructuredFiniteDifferenceDiscretizationCalculator().GetInterpolationWeightsVector(dpGrid, point, stencilSizes);
		}

		[[nodiscard]] Vector<CoordinateType> GetIntegrationWeightsVector(const Grid<Dimension, CoordinateType>& grid) const noexcept override
		{
			const auto& dpGrid = dynamic_cast<const DirectProductGrid<Dimension, CoordinateType>&>(grid);
			return StructuredFiniteDifferenceDiscretizationCalculator().GetIntegrationWeightsVector(dpGrid);
		}

	private:
		[[nodiscard]] static std::array<size_t, Dimension> FillStencils(size_t stencilSize) noexcept
		{
			std::array<size_t, Dimension> stencilSizes;
			std::fill(stencilSizes.begin(), stencilSizes.end(), stencilSize);
			return stencilSizes;
		}

		const std::array<size_t, Dimension> stencilSizes;
	};
}