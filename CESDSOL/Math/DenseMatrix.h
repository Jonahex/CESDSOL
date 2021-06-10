#pragma once

#include "Math/MultiLevelArray.h"

#include <cstddef>

namespace CESDSOL
{
	template<typename ScalarType>
	class DenseMatrix : private TwoLevelArray<ScalarType>
	{
	public:
		using TwoLevelArray<ScalarType>::begin;
		using TwoLevelArray<ScalarType>::end;
		using TwoLevelArray<ScalarType>::size;
		using TwoLevelArray<ScalarType>::data;
		using TwoLevelArray<ScalarType>::operator[];
		using TwoLevelArray<ScalarType>::Flatten;

		[[nodiscard]] constexpr size_t RowCount() const noexcept
		{
			return this->size();
		}

		[[nodiscard]] constexpr size_t ColumnCount() const noexcept
		{
			return this->operator[](0).size();
		}

		constexpr DenseMatrix(size_t rowCount, size_t columnCount):
			TwoLevelArray<ScalarType>({rowCount, columnCount})
		{}
	};
}