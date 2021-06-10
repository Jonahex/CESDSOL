#pragma once

#include "Math/LinearAlgebra.h"

namespace CESDSOL
{	
	template<Concepts::CSRMatrix MatrixType>
	[[nodiscard]] constexpr MatrixType MakeIdentityMatrix(size_t size)
	{
		auto result = MatrixType(size, size, size);
		for (size_t i = 0; i < size; i++)
		{
			result.SetRowCount(i, i);
			result.SetColumnIndex(i, i);
			result.SetValue(i, 1);
		}
		return result;
	}

	template<Concepts::CSRMatrix MatrixType>
	std::ostream& operator<<(std::ostream& stream, const MatrixType& matrix) noexcept
	{
		const auto rowCount = matrix.RowCount();
		const auto columnCount = matrix.ColumnCount();
		for (size_t i = 0; i < rowCount; i++)
		{
			size_t currentColumn = 0;
			for (size_t j = matrix.GetRowCount(i); j < matrix.GetRowCount(i + 1); ++j)
			{
				const auto columnIndex = matrix.GetColumnIndex(j);
				while (currentColumn < columnIndex)
				{
					stream << '0' << ' ';
					currentColumn++;
				}
				stream << matrix.GetValue(j) << ' ';
				currentColumn++;
			}
			while (currentColumn < columnCount)
			{
				stream << '0' << ' ';
				currentColumn++;
			}
			stream << '\n';
		}

		return stream;
	}
}