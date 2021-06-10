#pragma once

#include "Grid/Grid.h"
#include "Math/LinearAlgebra.h"
#include "Serialization/DataType.h"

namespace CESDSOL
{
	template<typename CoordinateType>
	struct SingleDimensionalGrid
	{
		Vector<CoordinateType> points;
		std::optional<CoordinateType> period;
	};

	template<size_t DimensionArg, typename CoordinateTypeArg>
	class DirectProductGrid : public Grid<DimensionArg, CoordinateTypeArg>
	{
	private:
		using Base = Grid<DimensionArg, CoordinateTypeArg>;

	public:
		using typename Base::CoordinateType;
		static constexpr size_t Dimension = Base::Dimension;

	private:
		std::array<SingleDimensionalGrid<CoordinateType>, Dimension> dimensionGrids;

		[[nodiscard]] static size_t GetDimensionSize(const std::array<SingleDimensionalGrid<CoordinateType>, Dimension>& dimensionGrids, size_t index) noexcept
		{
			return dimensionGrids[index].points.size();
		}

		static void GetMultiIndexBySingleIndex(const std::array<SingleDimensionalGrid<CoordinateType>, Dimension>& dimensionGrids, size_t size,
			size_t index, std::array<size_t, Dimension>& point) noexcept
		{
			if constexpr (Dimension == 1)
			{
				point[0] = index;
			}
			else if constexpr (Dimension == 2)
			{
				point[0] = index / GetDimensionSize(dimensionGrids, 1);
				point[1] = index % GetDimensionSize(dimensionGrids, 1);
			}
			else
			{
				size_t current = index;
				size_t tmp = size;
				for (size_t i = 0; i < Dimension; i++)
				{
					tmp /= GetDimensionSize(dimensionGrids, i);
					point[i] = current / tmp;
					current %= tmp;
				}
			}
		}

		[[nodiscard]] static Vector<Base::Point> MakePoints(const std::array<SingleDimensionalGrid<CoordinateType>, Dimension>& dimensionGrids) noexcept
		{
			const auto size = std::accumulate(dimensionGrids.begin(), dimensionGrids.end(), 1, 
				[](const auto product, const auto value) {return product * value.points.size(); });
			Vector<Base::Point> result(size);
			std::array<size_t, Dimension> multiIndex;
			for (size_t pointIndex = 0; pointIndex < size; ++pointIndex)
			{
				GetMultiIndexBySingleIndex(dimensionGrids, size, pointIndex, multiIndex);
				for (size_t dimensionIndex = 0; dimensionIndex < Dimension; ++dimensionIndex)
				{
					result[pointIndex].coordinates[dimensionIndex] = dimensionGrids[dimensionIndex].points[multiIndex[dimensionIndex]];
				}
				result[pointIndex].regionIndex = 0;
				for (int64_t dimensionIndex = Dimension - 1; dimensionIndex >= 0; --dimensionIndex)
				{
					if (!dimensionGrids[dimensionIndex].period.has_value())
					{
						if (multiIndex[dimensionIndex] == 0)
						{
							result[pointIndex].regionIndex = 2 * dimensionIndex + 1;
							break;
						}
						else if (multiIndex[dimensionIndex] == GetDimensionSize(dimensionGrids, dimensionIndex) - 1)
						{
							result[pointIndex].regionIndex = 2 * dimensionIndex + 2;
							break;
						}
					}
				}
			}
			return result;
		}

	public:		
		template<typename... GridRefTypes> requires((std::convertible_to<std::remove_cvref_t<GridRefTypes>, SingleDimensionalGrid<CoordinateType>> && ...) && sizeof...(GridRefTypes) == Dimension)
		DirectProductGrid(GridRefTypes&&... dimensionGrids) noexcept
			: Base(MakePoints(MakeArray(dimensionGrids...)), 1 + 2 * ((dimensionGrids.period.has_value() ? 0 : 1) + ...))
			, dimensionGrids{ std::forward<SingleDimensionalGrid<CoordinateType>>(dimensionGrids)... }
		{}

		template<typename GridRefType> requires(std::convertible_to<std::remove_cvref_t<GridRefType>, std::array<SingleDimensionalGrid<CoordinateType>, Dimension>>)
		DirectProductGrid(GridRefType&& dimensionGrids) noexcept
			: Base(MakePoints(dimensionGrids), std::reduce(dimensionGrids.begin(), dimensionGrids.end(), 1, [](const auto& sum, const auto& element) { return sum + 2 * !element.period.has_value(); }))
			, dimensionGrids{ std::forward<std::array<SingleDimensionalGrid<CoordinateType>, Dimension>>(dimensionGrids) }
		{}

		[[nodiscard]] const Vector<CoordinateType>& GetGrid(size_t index) const noexcept
		{
			return dimensionGrids[index].points;
		}

		[[nodiscard]] size_t GetDimensionSize(size_t index) const noexcept
		{
			return GetDimensionSize(dimensionGrids, index);
		}

		void GetMultiIndexBySingleIndex(size_t index, std::array<size_t, Dimension>& point) const noexcept
		{
			GetMultiIndexBySingleIndex(dimensionGrids, this->GetSize(), index, point);
		}

		[[nodiscard]] std::array<size_t, Dimension> GetMultiIndexBySingleIndex(size_t index) const noexcept
		{
			std::array<size_t, Dimension> result;
			GetMultiIndexBySingleIndex(index, result);
			return result;
		}

		[[nodiscard]] size_t GetSingleIndexByMultiIndex(const std::array<size_t, Dimension>& index) const noexcept
		{
			if constexpr (Dimension == 1)
			{
				return index[0];
			}
			else if constexpr (Dimension == 2)
			{
				return index[0] * GetDimensionSize(1) + index[1];
			}
			else if constexpr (Dimension == 3)
			{
				return (index[0] * GetDimensionSize(1) + index[1]) * GetDimensionSize(2) + index[2];
			}
			else
			{
				size_t result = index[Dimension - 1];
				size_t tmp = 1;
				for (i64 i = Dimension - 2; i >= 0; i--)
				{
					tmp *= GetDimensionSize(i + 1);
					result += index[i] * tmp;
				}
				return result;
			}
		}

		void GetCoordinatesByMultiIndex(const std::array<size_t, Dimension>& index, std::array<CoordinateType, Dimension>& point) const noexcept
		{
			for (size_t i = 0; i < Dimension; i++)
			{
				point[i] = dimensionGrids[i].points[index[i]];
			}
		}

		[[nodiscard]] std::array<CoordinateType, Dimension> GetCoordinatesByMultiIndex(const std::array<size_t, Dimension>& index) const noexcept
		{
			std::array<CoordinateType, Dimension> result;
			GetCoordinatesByMultiIndex(index, result);
			return result;
		}

		[[nodiscard]] bool IsPeriodicDimension(size_t index) const noexcept
		{
			return dimensionGrids[index].period.has_value();
		}

		[[nodiscard]] std::optional<CoordinateType> GetPeriod(size_t index) const noexcept
		{
			return dimensionGrids[index].period;
		}
		
	private:
		static constexpr Serializer::GridTypeData TypeData = { 1, static_cast<uint64_t>(Dimension), static_cast<uint64_t>(Serializer::ToDataType<CoordinateType>()) };

		[[nodiscard]] uint64_t GetGridType() const noexcept override
		{
			return TypeData.type;
		}

		struct Header
		{
			std::array<uint64_t, Dimension> dimensionSizes;
			std::array<CoordinateType, Dimension> periods;
		};

		[[nodiscard]] Array<uint8_t> GetHeader() const noexcept override
		{
			Header header;
			for (size_t dimensionIndex = 0; dimensionIndex < Dimension; ++dimensionIndex)
			{
				header.dimensionSizes[dimensionIndex] = GetDimensionSize(dimensionIndex);
				header.periods[dimensionIndex] = GetPeriod(dimensionIndex).has_value() ? *GetPeriod(dimensionIndex) : 0;
			}
			Array<uint8_t> result(sizeof(Header), reinterpret_cast<uint8_t*>(&header));
			return result;
		}

		[[nodiscard]] virtual Array<uint8_t> GetData() const noexcept
		{
			const auto pointCount =
				std::reduce(dimensionGrids.cbegin(), dimensionGrids.cend(), 0, [](auto sum, const auto& item) {return sum + item.points.size(); });
			Array<uint8_t> result = Array<uint8_t>(sizeof(CoordinateType) * pointCount);
			uint8_t* data = result.data();
			for (size_t dimensionIndex = 0; dimensionIndex < Dimension; ++dimensionIndex)
			{
				const auto diff = sizeof(CoordinateType) * GetDimensionSize(dimensionIndex);
				std::memcpy(data, dimensionGrids[dimensionIndex].points.data(), diff);
				data += diff;
			}
			return result;
		}

		[[nodiscard]] static std::shared_ptr<DirectProductGrid<Dimension, CoordinateType>>
			Load(const Array<uint8_t>& headerData, const Array<uint8_t>& data) noexcept
		{
			const auto header = reinterpret_cast<const Header*>(headerData.data());
			std::array<SingleDimensionalGrid<CoordinateType>, Dimension> grids;
			auto dataPtr = reinterpret_cast<const CoordinateType*>(data.data());
			for (size_t dimensionIndex = 0; dimensionIndex < Dimension; ++dimensionIndex)
			{
				grids[dimensionIndex].period = header->periods[dimensionIndex] == 0 ? std::optional<CoordinateType>() : std::optional<CoordinateType>(header->periods[dimensionIndex]);
				grids[dimensionIndex].points = Vector<CoordinateType>(header->dimensionSizes[dimensionIndex], dataPtr);
				dataPtr += header->dimensionSizes[dimensionIndex];
			}
			return std::make_shared<DirectProductGrid<Dimension, CoordinateType>>(std::move(grids));
		}

		RegisterGridLoader(TypeData, &Load);
	};
}