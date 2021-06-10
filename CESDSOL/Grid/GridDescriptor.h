#pragma once

#include <compare>

namespace CESDSOL
{
	template<size_t DimensionArg, typename CoordinateTypeArg = double>
	class GridDescriptor
	{
	public:
		static constexpr size_t Dimension = DimensionArg;
		using CoordinateType = CoordinateTypeArg;

		GridDescriptor(size_t aRegionCount) noexcept
			: regionCount(aRegionCount)
		{}

		[[nodiscard]] size_t GetRegionCount() const noexcept
		{
			return regionCount;
		}

		auto operator<=>(const GridDescriptor& second) const noexcept = default;

	private:
		size_t regionCount;
	};
}