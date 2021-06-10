#pragma once

#include "Concepts.h"
#include "Utils/Logger.h"
#include "Utils/Utils.h"

#include <concepts>

namespace CESDSOL
{
	template<Concepts::SparseVector SparseVectorType, Concepts::Vector VectorType> 
		requires std::common_with<typename SparseVectorType::value_type, typename VectorType::value_type>
	[[nodiscard]] constexpr auto DotProduct(const SparseVectorType& left, const VectorType& right)
	{
		typename VectorType::value_type result = 0;
		AssertE(left.ElementCount() == right.size(), MessageTag::Math, "Trying to calculate dot product for vectors of unequal size.");

		for (size_t i = 0; i < left.NonZeroCount(); i++)
		{
			result += left.GetValue(i) * right[left.GetIndex(i)];
		}
		return result;
	}
		
	template<Concepts::SparseVector SparseVectorType, Concepts::Vector VectorType> 
		requires std::common_with<typename SparseVectorType::value_type, typename VectorType::value_type>
	[[nodiscard]] constexpr auto DotProduct(const VectorType& left, const SparseVectorType& right)
	{
		return DotProduct(right, left);
	}

	template<Concepts::SparseVector SparseVectorType>
	[[nodiscard]] constexpr SparseVectorType DirectProductAsVector(const SparseVectorType& left, const SparseVectorType& right)
	{
		auto result = SparseVectorType(left.ElementCount() * right.ElementCount(), left.NonZeroCount() * right.NonZeroCount());
		size_t k = 0;
		for (size_t i = 0; i < left.NonZeroCount(); i++)
		{
			for (size_t j = 0; j < right.NonZeroCount(); j++)
			{
				result.SetValue(k, left.GetValue(i) * right.GetValue(j));
				result.SetIndex(k++, left.GetIndex(i) * right.ElementCount() + right.GetIndex(j));
			}
		}
		return result;
	}
}