#pragma once

#include "Grid/BaseGrid.h"
#include "Grid/GridDescriptor.h"
#include "Math/LinearAlgebra.h"
#include "Serialization/DataType.h"
#include "Serialization/GridCache.h"

namespace CESDSOL
{
	template<size_t DimensionArg, typename CoordinateTypeArg>
	class Grid : public BaseGrid
	{
	public:
		using CoordinateType = CoordinateTypeArg;
		static constexpr size_t Dimension = DimensionArg;

		struct Point
		{
			std::array<CoordinateType, Dimension> coordinates;
			size_t regionIndex;
		};

		struct SerializedData
		{
			Serializer::GridTypeData type;
			Array<uint8_t> header;
			Array<uint8_t> data;
		};

		[[nodiscard]] size_t GetSize() const noexcept
		{
			return points.size();
		}

		[[nodiscard]] size_t GetRegionIndex(size_t pointIndex) const noexcept
		{
			return points[pointIndex].regionIndex;
		}

		[[nodiscard]] const std::array<CoordinateType, Dimension>& GetCoordinates(size_t pointIndex) const noexcept
		{
			return points[pointIndex].coordinates;
		}

		[[nodiscard]] size_t GetRegionCount() const noexcept
		{
			return regionCount;
		}

		[[nodiscard]] GridDescriptor<Dimension, CoordinateType> GetDescriptor() const noexcept
		{
			return { GetRegionCount() };
		}

		[[nodiscard]] SerializedData Save() const noexcept
		{
			SerializedData result;
			result.type = TypeData;
			result.type.type = GetGridType();
			result.header = GetHeader();
			result.data = GetData();
			return result;
		}

		template<typename VectorRefType> requires (std::same_as<std::remove_cvref_t<VectorRefType>, Vector<Point>>)
		Grid(VectorRefType&& aPoints, size_t aRegionCount) noexcept
			: points(std::forward<Vector<Point>>(aPoints))
			, regionCount(aRegionCount)
		{}

		virtual ~Grid() = default;

	private:
		struct Header
		{
			size_t regionCount;
			size_t pointCount;
		};

	protected:
		[[nodiscard]] virtual uint64_t GetGridType() const noexcept
		{
			return TypeData.type;
		}

		[[nodiscard]] virtual Array<uint8_t> GetHeader() const noexcept
		{
			Header header{ regionCount, GetSize() };
			Array<uint8_t> result(sizeof(Header), reinterpret_cast<uint8_t*>(&header));
			return result;
		}

		[[nodiscard]] virtual Array<uint8_t> GetData() const noexcept
		{
			Array<uint8_t> result = Array<uint8_t>(sizeof(Point) * points.size(), reinterpret_cast<uint8_t*>(points.data()));
			return result;
		}

	private:
		static constexpr Serializer::GridTypeData TypeData = { 0, static_cast<uint64_t>(Dimension), static_cast<uint64_t>(Serializer::ToDataType<CoordinateType>()) };

		[[nodiscard]] static std::shared_ptr<Grid<Dimension, CoordinateType>> Load(const Array<uint8_t>& headerData, const Array<uint8_t>& data) noexcept
		{
			const auto header = reinterpret_cast<const Header*>(headerData.data());
			Vector<Point> points(header->pointCount, reinterpret_cast<const Point*>(data.data()));
			return std::make_shared<Grid<Dimension, CoordinateType>>(points, header->regionCount);
		}

		Vector<Point> points;
		size_t regionCount;

		RegisterGridLoader(TypeData, &Load);
	};

	template<typename GridType, typename FunctorType>
	[[nodiscard]] auto Apply(const GridType& grid, const FunctorType& functor) noexcept
	{
		using ScalarType = typename GridType::CoordinateType;
		auto result = Vector<ScalarType>(grid.GetSize());
		for (size_t index = 0; index < grid.GetSize(); ++index)
		{
			result[index] = functor(grid.GetCoordinates(index));
		}
		return result;
	}

	template<size_t Dimension, typename CoordinateType>
	[[nodiscard]] bool operator==(const Grid<Dimension, CoordinateType>& first, const Grid<Dimension, CoordinateType>& second) noexcept
	{
		if (first.GetSize() != second.GetSize())
		{
			return false;
		}
		for (size_t pointIndex = 0; pointIndex < first.GetSize(); ++pointIndex)
		{
			if (first.GetCoordinates(pointIndex) != second.GetCoordinates(pointIndex))
			{
				return false;
			}
		}
		return true;
	}
}