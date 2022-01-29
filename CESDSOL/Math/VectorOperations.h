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

	template<Concepts::Vector VectorType>
	[[nodiscard]] constexpr auto Transform(const VectorType& vector, const std::function<typename VectorType::value_type(typename VectorType::value_type)>& functor) noexcept
	{
		auto result = Vector<typename VectorType::value_type>(vector.size());
		std::transform(vector.begin(), vector.end(), result.begin(), functor);
		return result;
	}

	template<Concepts::Vector VectorType, typename ScalarType>
	[[nodiscard]] constexpr auto LowerBoundIndexBinary(const VectorType& array, const ScalarType& value) noexcept
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

	template<typename ScalarType, Concepts::Vector XVectorType, Concepts::Vector YVectorType>
	constexpr void AXPY(ScalarType a, const XVectorType& x, YVectorType& y) noexcept
	{
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to perform AXPY operation with vectors of unequal size.");
		LinearAlgebra::AXPY(a, x.data(), y.data(), x.size());
	}

	template<typename AScalarType, Concepts::Vector XVectorType, typename BScalarType, Concepts::Vector YVectorType>
	constexpr void AXPBY(typename AScalarType a, const XVectorType& x, BScalarType b, YVectorType& y) noexcept
	{
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to perform AXPBY operation with vectors of unequal size.");
		LinearAlgebra::AXPBY(a, x.data(), b, y.data(), x.size());
	}

	template<typename ScalarType, Concepts::Vector VectorType>
	constexpr void Scale(typename ScalarType a, VectorType& x) noexcept
	{
		LinearAlgebra::Scale(a, x.data(), x.size());
	}

	template<Concepts::Vector LeftVectorType, Concepts::Vector RightVectorType>
	[[nodiscard]] constexpr auto DirectProductAsVector(const LeftVectorType& left, const RightVectorType& right) noexcept
	{
		using ScalarType = std::common_type_t<typename LeftVectorType::value_type, typename RightVectorType::value_type>;
		auto result = Vector<ScalarType>(left.size() * right.size());
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

	template<Concepts::Vector VectorType>
	constexpr auto Norm2(const VectorType& x) noexcept
	{
		return LinearAlgebra::Norm2(x.data(), x.size());
	}

	template<Concepts::Vector LeftVectorType, Concepts::Vector RightVectorType>
	[[nodiscard]] auto Add(const LeftVectorType& left, const RightVectorType& right) noexcept
	{
		AssertE(left.size() == right.size(), MessageTag::Math, "Trying to add vectors of unequal size.");
		using ScalarType = std::common_type_t<typename LeftVectorType::value_type, typename RightVectorType::value_type>;
		auto result = Vector<ScalarType>(left.size());
		LinearAlgebra::Add(left.data(), right.data(), result.data(), left.size());
		return result;
	}

	template<Concepts::Vector LeftVectorType, Concepts::Vector RightVectorType>
	[[nodiscard]] auto Subtract(const LeftVectorType& left, const RightVectorType& right) noexcept
	{
		AssertE(left.size() == right.size(), MessageTag::Math, "Trying to subtract vectors of unequal size.");
		using ScalarType = std::common_type_t<typename LeftVectorType::value_type, typename RightVectorType::value_type>;
		auto result = Vector<ScalarType>(left.size());
		LinearAlgebra::Subtract(left.data(), right.data(), result.data(), left.size());
		return result;
	}

	template<Concepts::Vector LeftVectorType, Concepts::Vector RightVectorType>
	[[nodiscard]] auto Multiply(const LeftVectorType& left, const RightVectorType& right) noexcept
	{
		AssertE(left.size() == right.size(), MessageTag::Math, "Trying to multiply vectors of unequal size.");
		using ScalarType = std::common_type_t<typename LeftVectorType::value_type, typename RightVectorType::value_type>;
		auto result = Vector<ScalarType>(left.size());
		LinearAlgebra::Multiply(left.data(), right.data(), result.data(), left.size());
		return result;
	}

	template<Concepts::Vector AVectorType, Concepts::Vector BVectorType>
	[[nodiscard]] auto DotProduct(const AVectorType& a, const BVectorType& b) noexcept
	{
		AssertE(a.size() == b.size(), MessageTag::Math, "Trying to find dot product of vectors of unequal size.");
		return LinearAlgebra::DotProduct(a.data(), b.data(), a.size());
	}

	template<Concepts::Vector VectorType, typename ScalarType>
	void Fill(VectorType& vector, ScalarType value) noexcept
	{
		return LinearAlgebra::Fill(vector.data(), vector.size(), value);
	}
}