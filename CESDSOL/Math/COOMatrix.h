#pragma once

#include "Array.h"

#include <cstddef>
#include <vector>

namespace CESDSOL::Native
{
	template<typename ScalarType, typename IndexType = size_t, IndexType StartingIndex = 0>
	class COOMatrix
	{
	private:
		Array<ScalarType> values;
		Array<IndexType> rows;
		Array<IndexType> columns;

		size_t rowCount = 0;
		size_t columnCount = 0;

	public:
		using value_type = ScalarType;
		using size_type = IndexType;
		using IndexType = IndexType;
		static constexpr IndexType StartingIndex = StartingIndex;

		[[nodiscard]] constexpr size_t RowCount() const noexcept
		{
			return rowCount;
		}

		[[nodiscard]] constexpr size_t ColumnCount() const noexcept
		{
			return columnCount;
		}

		[[nodiscard]] constexpr size_t ElementCount() const noexcept
		{
			return rowCount * columnCount;
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

		[[nodiscard]] constexpr const Array<IndexType>& GetColumns() const noexcept
		{
			return columns;
		}

		[[nodiscard]] constexpr Array<IndexType>& GetColumns() noexcept
		{
			return columns;
		}

		[[nodiscard]] constexpr const Array<IndexType>& GetRows() const noexcept
		{
			return rows;
		}

		[[nodiscard]] constexpr Array<IndexType>& GetRows() noexcept
		{
			return rows;
		}

		constexpr void SetValue(size_t index, ScalarType value) noexcept
		{
			values[index] = value;
		}

		[[nodiscard]] constexpr ScalarType GetValue(size_t index) const noexcept
		{
			return values[index];
		}

		constexpr void SetColumn(size_t index, IndexType value) noexcept
		{
			columns[index] = value + StartingIndex;
		}

		[[nodiscard]] constexpr IndexType GetColumn(size_t index) const noexcept
		{
			return columns[index] - StartingIndex;
		}

		constexpr void SetRow(size_t index, IndexType value) noexcept
		{
			rows[index] = value + StartingIndex;
		}

		[[nodiscard]] constexpr IndexType GetRow(size_t index) const noexcept
		{
			return rows[index] - StartingIndex;
		}

		constexpr COOMatrix() noexcept = default;

		constexpr COOMatrix(size_t rowCount, size_t columnCount, size_t nonZeroCount) noexcept :
			values(nonZeroCount),
			columns(nonZeroCount),
			rows(nonZeroCount),
			rowCount(rowCount),
			columnCount(columnCount)
		{}
	};
}