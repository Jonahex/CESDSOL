#pragma once

#include "Math/Concepts.h"
#include "Math/LinearAlgebra.h"

#include <algorithm>
#include <concepts>

namespace CESDSOL
{
	template<typename ScalarType, Concepts::Vector ArrayType = Vector<ScalarType>> 
		requires std::convertible_to<ScalarType, typename ArrayType::value_type>
	[[nodiscard]] constexpr ArrayType MakeUniformRange(ScalarType first, ScalarType last, size_t size) noexcept
	{
		auto range = ArrayType(size + 1);
		ScalarType current = first;
		ScalarType step = (last - first) / size;
		for (auto& x : range)
		{
			x = current;
			current += step;
		}
		return range;
	}

	template<typename ScalarType> 
	[[nodiscard]] constexpr Vector<ScalarType> Transform(const Vector<ScalarType>& vector, const std::function<ScalarType(ScalarType)>& functor) noexcept
	{
		auto result = Vector<ScalarType>(vector.size());
		std::transform(vector.begin(), vector.end(), result.begin(), functor);
		return result;
	}

	template<Concepts::Vector ArrayType, typename ScalarType>
	[[nodiscard]] constexpr auto LowerBoundIndexBinary(const ArrayType& array, const ScalarType& value) noexcept
	{
		const auto lowerBound = std::lower_bound(array.begin(), array.end(), value);
		return (lowerBound - array.begin());
	}

	template<Concepts::Vector SourceArrayType, Concepts::Vector DestinationArrayType>
	constexpr void Copy(const SourceArrayType& x, DestinationArrayType& y) noexcept
	{
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to copy vectors with unequal size.");
		LinearAlgebra::Copy(x.data(), y.data(), x.size());
	}

	template<Concepts::Vector ArrayType1, Concepts::Vector ArrayType2>
	constexpr void AXPY(typename ArrayType1::value_type a, const ArrayType1& x, ArrayType2& y) noexcept
	{
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to perform AXPY operation with vectors of unequal size.");
		LinearAlgebra::AXPY(a, x.data(), y.data(), x.size());
	}

	template<Concepts::Vector ArrayType1, Concepts::Vector ArrayType2>
	constexpr void AXPBY(typename ArrayType1::value_type a, const ArrayType1& x, typename ArrayType1::value_type b, ArrayType2& y) noexcept
	{
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to perform AXPBY operation with vectors of unequal size.");
		LinearAlgebra::AXPBY(a, x.data(), b, y.data(), x.size());
	}

	template<Concepts::Vector ArrayType>
	constexpr void Scale(typename ArrayType::value_type a, ArrayType& x) noexcept
	{
		LinearAlgebra::Scale(a, x.data(), x.size());
	}

	template<Concepts::Vector VectorType>
	[[nodiscard]] constexpr VectorType DirectProductAsVector(const VectorType& left, const VectorType& right) noexcept
	{
		auto result = VectorType(left.size() * right.size());
		size_t k = 0;
		for (size_t i = 0; i < left.size(); i++)
		{
			for (size_t j = 0; j < right.size(); j++)
			{
				result[i * right.size() + j] = left[i] * right[j];
			}
		}
		return result;
	}

	template<Concepts::Vector ArrayType>
	constexpr auto Norm2(const ArrayType& x) noexcept
	{
		return LinearAlgebra::Norm2(x.data(), x.size());
	}
}