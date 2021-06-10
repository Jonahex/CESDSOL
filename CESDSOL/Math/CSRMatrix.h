#pragma once

#include "Array.h"

#include <cstddef>

namespace CESDSOL::Native
{
	template<typename ScalarType, typename IndexType = size_t, IndexType StartingIndexValue = 1>
	class CSRMatrix
	{
	protected:
		Array<ScalarType> values;
		Array<IndexType> columnIndices;
		Array<IndexType> rowCounts;

		size_t rowCount = 0;
		size_t columnCount = 0;

	public:
		using value_type = ScalarType;
		using size_type = size_t;
		using index_type = IndexType;
		static constexpr IndexType StartingIndex = StartingIndexValue;

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
			return rowCounts[rowCount] - StartingIndex;
		}

		[[nodiscard]] constexpr const Array<ScalarType>& GetValues() const noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr Array<ScalarType>& GetValues() noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr const Array<IndexType>& GetColumnIndices() const noexcept
		{
			return columnIndices;
		}

		[[nodiscard]] constexpr Array<IndexType>& GetColumnIndices() noexcept
		{
			return columnIndices;
		}

		[[nodiscard]] constexpr const Array<IndexType>& GetRowCounts() const noexcept
		{
			return rowCounts;
		}

		[[nodiscard]] constexpr Array<IndexType>& GetRowCounts() noexcept
		{
			return rowCounts;
		}

		constexpr void SetValue(size_t index, ScalarType value) noexcept
		{
			values[index] = value;
		}

		[[nodiscard]] constexpr ScalarType GetValue(size_t index) const noexcept
		{
			return values[index];
		}

		constexpr void SetColumnIndex(size_t index, IndexType value) noexcept
		{
			columnIndices[index] = value + StartingIndex;
		}

		[[nodiscard]] constexpr IndexType GetColumnIndex(size_t index) const noexcept
		{
			return columnIndices[index] - StartingIndex;
		}

		constexpr void SetRowCount(size_t index, IndexType value) noexcept
		{
			rowCounts[index] = value + StartingIndex;
		}

		[[nodiscard]] constexpr IndexType GetRowCount(size_t index) const noexcept
		{
			return rowCounts[index] - StartingIndex;
		}

		[[nodiscard]] constexpr IndexType GetRowLength(size_t index) const noexcept
		{
			return rowCounts[index + 1] - rowCounts[index];
		}

		constexpr void Nullify() noexcept
		{
			for (auto& value : values)
			{
				value = 0;
			}
		}

		void ReplaceValues(const Array<ScalarType>& newValues) noexcept
		{
			AssertE(this->values.size() == newValues.size(), MessageTag::Math,
				"Incompatible number of nonzero values in CSR matrix!");
			values = newValues;
		}

		constexpr CSRMatrix() noexcept = default;

		constexpr CSRMatrix(size_t aRowCount, size_t aColumnCount, size_t nonZeroCount) noexcept :
			values(nonZeroCount),
			columnIndices(nonZeroCount),
			rowCounts(aRowCount + 1),
			rowCount(aRowCount),
			columnCount(aColumnCount)
		{
			rowCounts[rowCount] = nonZeroCount + StartingIndex;
		}
	};
}