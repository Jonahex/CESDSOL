#pragma once

#include <cstddef>

namespace CESDSOL::Native
{
	template<typename ScalarType>
	constexpr void Copy(const ScalarType* source, ScalarType* destination, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			destination[i] = source[i];
		}
	}

	template<typename ScalarType>
	constexpr void AXPY(ScalarType a, const ScalarType* x, ScalarType* y, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			y[i] += a * x[i];
		}
	}

	template<typename ScalarType>
	constexpr void AXPBY(ScalarType a, const ScalarType* x, ScalarType b, ScalarType* y, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			y[i] = a * x[i] + b * y[i];
		}
	}

	template<typename ScalarType>
	constexpr void Scale(ScalarType a, ScalarType* x, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			x[i] *= a;
		}
	}

	template<typename ScalarType>
	constexpr ScalarType Norm2(const ScalarType* x, size_t count) noexcept
	{
		ScalarType result = 0;
		for (size_t i = 0; i < count; i++)
		{
			result += x[i] * x[i];
		}
		return std::sqrt(result);
	}

	template<typename ScalarType>
	constexpr ScalarType Norm2(const std::complex<ScalarType>* x, size_t count) noexcept
	{
		ScalarType result = 0;
		for (size_t i = 0; i < count; i++)
		{
			result += x[i] * std::conj(x[i]);
		}
		return std::sqrt(result);
	}
}