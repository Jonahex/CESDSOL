#pragma once

#include "Grid/Grid.h"

namespace CESDSOL
{
	template<size_t DimensionArg, template<typename> typename MatrixTypeArg = CSRMatrix, typename CoordinateTypeArg = double>
	class Discretization
	{
	public:
		static constexpr size_t Dimension = DimensionArg;
		using CoordinateType = CoordinateTypeArg;
		using MatrixType = MatrixTypeArg<CoordinateType>;

		[[nodiscard]] virtual MatrixType GetDifferentiationMatrix(const Grid<Dimension, CoordinateType>& grid,
			const std::array<size_t, Dimension>& derivativeOrders) const noexcept = 0;
		[[nodiscard]] virtual SparseVector<CoordinateType> GetInterpolationWeightsVector(const Grid<Dimension, CoordinateType>& grid,
			const std::array<CoordinateType, Dimension>& point) const noexcept = 0;
		[[nodiscard]] virtual Vector<CoordinateType> GetIntegrationWeightsVector(const Grid<Dimension, CoordinateType>& grid) const noexcept = 0;

		virtual ~Discretization() = default;
	};
}