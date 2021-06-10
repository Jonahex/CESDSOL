#pragma once

#include "Array.h"

namespace CESDSOL
{
	template<typename ScalarType, typename IndexType = size_t>
	class SparseVector
	{
	private:
		Array<ScalarType> values;
		Array<IndexType> indices;

		size_t length = 0;

	public:
		using value_type = ScalarType;
		using size_type = size_t;
		using index_type = IndexType;

		[[nodiscard]] constexpr size_t ElementCount() const noexcept
		{
			return length;
		}

		[[nodiscard]] constexpr size_t NonZeroCount() const noexcept
		{
			return values.size();
		}

		[[nodiscard]] constexpr const Array<ScalarType>& GetValues() const noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr Array<ScalarType>& GetValues() noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr const Array<IndexType>& GetIndices() const noexcept
		{
			return indices;
		}

		[[nodiscard]] constexpr Array<IndexType>& GetIndices() noexcept
		{
			return indices;
		}

		constexpr void SetValue(size_t index, ScalarType value) noexcept
		{
			values[index] = value;
		}

		[[nodiscard]] constexpr ScalarType GetValue(size_t index) const noexcept
		{
			return values[index];
		}

		constexpr void SetIndex(size_t index, IndexType value) noexcept
		{
			indices[index] = value;
		}

		[[nodiscard]] constexpr IndexType GetIndex(size_t index) const noexcept
		{
			return indices[index];
		}

		constexpr SparseVector() noexcept = default;

		constexpr SparseVector(size_t length, size_t nonZeroCount) noexcept:
			values(nonZeroCount),
			indices(nonZeroCount),
			length(length)
		{}
	};
}