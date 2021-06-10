#pragma once

#include <concepts>
#include <cstddef>

namespace CESDSOL::Concepts
{
	template<typename ArrayType>
	concept Vector = requires(ArrayType array)
	{
		typename ArrayType::value_type;
		typename ArrayType::size_type;

		{array.data()} -> std::convertible_to<typename ArrayType::value_type*>;
		{array.size()} -> std::convertible_to<typename ArrayType::size_type>;
	};

	template<typename ArrayType>
	concept SparseVector = requires(ArrayType array, size_t index,
		typename ArrayType::value_type value, typename ArrayType::index_type inputIndex)
	{
		typename ArrayType::value_type;
		typename ArrayType::size_type;
		typename ArrayType::index_type;

		{array.ElementCount()} -> std::same_as<size_t>;
		{array.NonZeroCount()} -> std::same_as<size_t>;
		{array.GetValue(index)} -> std::same_as<typename ArrayType::value_type>;
		{array.GetIndex(index)} -> std::same_as<typename ArrayType::index_type>;
		{array.SetValue(index, value)} -> std::same_as<void>;
		{array.SetIndex(index, inputIndex)} -> std::same_as<void>;
	};

	template<typename MatrixType>
	concept CSRMatrix = requires(MatrixType matrix, size_t index, 
		typename MatrixType::value_type value, typename MatrixType::index_type inputIndex)
	{
		typename MatrixType::value_type;
		typename MatrixType::size_type;
		typename MatrixType::index_type;

		{matrix.RowCount()} -> std::same_as<size_t>;
		{matrix.ColumnCount()} -> std::same_as<size_t>;
		{matrix.NonZeroCount()} -> std::same_as<size_t>;
		{matrix.GetValue(index)} -> std::same_as<typename MatrixType::value_type>;
		{matrix.GetColumnIndex(index)} -> std::same_as<typename MatrixType::index_type>;
		{matrix.GetRowCount(index)} -> std::same_as<typename MatrixType::index_type>;
		{matrix.SetValue(index, value)} -> std::same_as<void>;
		{matrix.SetColumnIndex(index, inputIndex)} -> std::same_as<void>;
		{matrix.SetRowCount(index, inputIndex)} -> std::same_as<void>;
	};
}