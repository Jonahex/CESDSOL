#pragma once

#include <cstddef>

namespace CESDSOL::Native
{
	template<typename SourceScalarType, typename DestinationScalarType>
	constexpr void Copy(const SourceScalarType* source, DestinationScalarType* destination, size_t count) noexcept
	{
		std::copy(source, source + count, destination);
	}

	template<typename AScalarType, typename XScalarType, typename YScalarType>
	constexpr void AXPY(AScalarType a, const XScalarType* x, YScalarType* y, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			y[i] += a * x[i];
		}
	}

	template<typename AScalarType, typename XScalarType, typename BScalarType, typename YScalarType>
	constexpr void AXPBY(AScalarType a, const XScalarType* x, BScalarType b, YScalarType* y, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			y[i] = a * x[i] + b * y[i];
		}
	}

	template<typename AScalarType, typename XScalarType>
	constexpr void Scale(AScalarType a, XScalarType* x, size_t count) noexcept
	{
		for (size_t i = 0; i < count; i++)
		{
			x[i] *= a;
		}
	}

	template<typename ScalarType>
	[[nodiscard]] ScalarType Norm2(const ScalarType* x, size_t count) noexcept
	{
		ScalarType result = 0;
		for (size_t i = 0; i < count; i++)
		{
			result += x[i] * x[i];
		}
		return std::sqrt(result);
	}

	template<typename ScalarType>
	[[nodiscard]] ScalarType Norm2(const std::complex<ScalarType>* x, size_t count) noexcept
	{
		ScalarType result = 0;
		for (size_t i = 0; i < count; i++)
		{
			result += x[i] * std::conj(x[i]);
		}
		return std::sqrt(result);
	}

	template<typename AScalarType, typename BScalarType>
	void Add(const AScalarType* a, const BScalarType* b, std::common_type_t<AScalarType, BScalarType>* x, size_t count) noexcept
	{
		for (size_t i = 0; i < count; ++i)
		{
			x[i] = a[i] + b[i];
		}
	}

	template<typename AScalarType, typename BScalarType>
	void Subtract(const AScalarType* a, const BScalarType* b, std::common_type_t<AScalarType, BScalarType>* x, size_t count) noexcept
	{
		for (size_t i = 0; i < count; ++i)
		{
			x[i] = a[i] - b[i];
		}
	}

	template<typename AScalarType, typename BScalarType>
	void Multiply(const AScalarType* a, const BScalarType* b, std::common_type_t<AScalarType, BScalarType>* x, size_t count) noexcept
	{
		for (size_t i = 0; i < count; ++i)
		{
			x[i] = a[i] * b[i];
		}
	}

	template<typename AScalarType, typename BScalarType>
	[[nodiscard]] auto DotProduct(const AScalarType* a, const BScalarType* b, size_t count) noexcept
	{
		decltype(*a * *b) result = 0;
		for (size_t i = 0; i < count; ++i)
		{
			result += a[i] * b[i];
		}
		return result;
	}

	template<typename SourceScalarType, typename DestinationScalarType>
	void Fill(DestinationScalarType* destination, size_t count, SourceScalarType value) noexcept
	{
		std::fill(destination, destination + count, value);
	}
}