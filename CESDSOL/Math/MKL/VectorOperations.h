#pragma once

#include "Math/MKL/Concepts.h"
#include "Math/MKL/Utils.h"
#include "Math/Native/VectorOperations.h"

#include "mkl.h"

#include <cstddef>

namespace CESDSOL::MKL
{	
	template<MKLScalar ScalarType>
	constexpr void Copy(const ScalarType* x, ScalarType* y, size_t count) noexcept
	{
#define COPY_OPERATION(type) cblas_##type##copy
		MKL_CBLAS_OPERATION(COPY_OPERATION, static_cast<int>(count), x, 1, y, 1);
#undef COPY_OPERATION
	}
	using Native::Copy;

	template<MKLScalar ScalarType>
	constexpr void AXPY(ScalarType a, const ScalarType* x, ScalarType* y, size_t count) noexcept
	{
#define AXPY_OPERATION(type) cblas_##type##axpy
		MKL_CBLAS_OPERATION(AXPY_OPERATION, static_cast<int>(count), a, x, 1, y, 1);
#undef AXPY_OPERATION
	}
	using Native::AXPY;

	template<MKLScalar ScalarType>
	constexpr void AXPBY(ScalarType a, const ScalarType* x, ScalarType b, ScalarType* y, size_t count) noexcept
	{
#define AXPBY_OPERATION(type) cblas_##type##axpby
		MKL_CBLAS_OPERATION(AXPBY_OPERATION, static_cast<int>(count), a, x, 1, b, y, 1);
#undef AXPY_OPERATION
	}
	using Native::AXPBY;

	template<MKLScalar ScalarType>
	constexpr void Scale(ScalarType a, ScalarType* x, size_t count) noexcept
	{
#define AXPY_OPERATION(type) cblas_##type##scal
		MKL_CBLAS_OPERATION(AXPY_OPERATION, static_cast<int>(count), a, x, 1);
#undef AXPY_OPERATION
	}
	using Native::AXPY;

	template<MKLScalar ScalarType>
	constexpr ScalarType Norm2(const ScalarType* x, size_t count) noexcept
	{
		if constexpr (std::is_same_v<ScalarType, f32>)
		{
			return cblas_snrm2(static_cast<int>(count), x, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, f64>)
		{
			return cblas_dnrm2(static_cast<int>(count), x, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, c32>)
		{
			return cblas_scnrm2(static_cast<int>(count), static_cast<const void*>(x), 1);
		}
		else if constexpr (std::is_same_v<ScalarType, c64>)
		{
			return cblas_dznrm2(static_cast<int>(count), static_cast<const void*>(x), 1);
		}
	}
	using Native::Norm2;
}