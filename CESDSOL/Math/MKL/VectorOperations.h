#pragma once

#include "Math/MKL/Concepts.h"
#include "Math/MKL/Utils.h"
#include "Math/Native/VectorOperations.h"

#include "mkl.h"

#include <cstddef>

namespace CESDSOL::MKL
{	
	template<MKLScalar ScalarType>
	void Copy(const ScalarType* x, ScalarType* y, size_t count) noexcept
	{
#define COPY_OPERATION(type) cblas_##type##copy
		MKL_CBLAS_OPERATION(COPY_OPERATION, static_cast<MKL_INT>(count), x, 1, y, 1);
#undef COPY_OPERATION
	}
	using Native::Copy;

	template<MKLScalar ScalarType>
	void AXPY(ScalarType a, const ScalarType* x, ScalarType* y, size_t count) noexcept
	{
#define AXPY_OPERATION(type) cblas_##type##axpy
		MKL_CBLAS_OPERATION(AXPY_OPERATION, static_cast<MKL_INT>(count), a, x, 1, y, 1);
#undef AXPY_OPERATION
	}
	using Native::AXPY;

	template<MKLScalar ScalarType>
	void AXPBY(ScalarType a, const ScalarType* x, ScalarType b, ScalarType* y, size_t count) noexcept
	{
#define AXPBY_OPERATION(type) cblas_##type##axpby
		MKL_CBLAS_OPERATION(AXPBY_OPERATION, static_cast<MKL_INT>(count), a, x, 1, b, y, 1);
#undef AXPY_OPERATION
	}
	using Native::AXPBY;

	template<MKLScalar ScalarType>
	void Scale(ScalarType a, ScalarType* x, size_t count) noexcept
	{
#define AXPY_OPERATION(type) cblas_##type##scal
		MKL_CBLAS_OPERATION(AXPY_OPERATION, static_cast<MKL_INT>(count), a, x, 1);
#undef AXPY_OPERATION
	}
	using Native::Scale;

	template<MKLScalar ScalarType>
	[[nodiscard]] ScalarType Norm2(const ScalarType* x, size_t count) noexcept
	{
		if constexpr (std::is_same_v<ScalarType, f32>)
		{
			return cblas_snrm2(static_cast<MKL_INT>(count), x, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, f64>)
		{
			return cblas_dnrm2(static_cast<MKL_INT>(count), x, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, c32>)
		{
			return cblas_scnrm2(static_cast<MKL_INT>(count), static_cast<const void*>(x), 1);
		}
		else if constexpr (std::is_same_v<ScalarType, c64>)
		{
			return cblas_dznrm2(static_cast<MKL_INT>(count), static_cast<const void*>(x), 1);
		}
	}
	using Native::Norm2;

	template<MKLScalar ScalarType>
	[[nodiscard]] auto DotProduct(const ScalarType* a, const ScalarType* b, size_t count) noexcept
	{
		if constexpr (std::is_same_v<ScalarType, f32>)
		{
			return cblas_sdot(static_cast<MKL_INT>(count), a, 1, b, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, f64>)
		{
			return cblas_ddot(static_cast<MKL_INT>(count), a, 1, b, 1);
		}
		else if constexpr (std::is_same_v<ScalarType, c32>)
		{
			c32 result;
			cblas_cdotu_sub(static_cast<MKL_INT>(count), a, 1, b, 1, &result);
			return result;
		}
		else if constexpr (std::is_same_v<ScalarType, c64>)
		{
			c64 result;
			cblas_ddotu_sub(static_cast<MKL_INT>(count), a, 1, b, 1, &result);
			return result;
		}
	}
	using Native::DotProduct;

	template<typename ScalarType>
	void Add(const ScalarType* a, const ScalarType* b, ScalarType* x, size_t count) noexcept
	{
#define ADD_OPERATION(type) v##type##Add
		MKL_VML_OPERATION(ADD_OPERATION, static_cast<MKL_INT>(count), a, b, x);
#undef ADD_OPERATION
	}
	using Native::Add;

	template<typename ScalarType>
	void Subtract(const ScalarType* a, const ScalarType* b, ScalarType* x, size_t count) noexcept
	{
#define SUB_OPERATION(type) v##type##Sub
		MKL_VML_OPERATION(SUB_OPERATION, static_cast<MKL_INT>(count), a, b, x);
#undef SUB_OPERATION
	}
	using Native::Subtract;

	template<typename ScalarType>
	void Multiply(const ScalarType* a, const ScalarType* b, ScalarType* x, size_t count) noexcept
	{
#define MUL_OPERATION(type) v##type##Sub
		MKL_VML_OPERATION(MUL_OPERATION, static_cast<MKL_INT>(count), a, b, x);
#undef MUL_OPERATION
	}
	using Native::Multiply;

	using Native::Fill;
}