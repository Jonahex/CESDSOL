#pragma once

namespace CESDSOL::MKL
{
	static constexpr const char* SparseStatusMessages[]
	{
		"the operation was successful.",
		"the routine encountered an empty handle or matrix array.",
		"internal memory allocation failed.",
		"the input parameters contain an invalid value.",
		"execution failed.",
		"an error in algorithm implementation occurred.",
		"the requested operation is not supported."
	};
}

#define MKL_CBLAS_OPERATION(operation, ...)														\
	if constexpr (std::same_as<ScalarType, f32>)												\
	{																							\
		operation(s)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, f64>)											\
	{																							\
		operation(d)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, c32>)											\
	{																							\
		operation(c)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, c64>)											\
	{																							\
		operation(z)(__VA_ARGS__);																\
	}

#define MKL_SPARSE_OPERATION(operation, name, ...)												\
	sparse_status_t status;																		\
	if constexpr (std::same_as<ScalarType, f32>)												\
	{																							\
		status = operation(s)(__VA_ARGS__);														\
	}																							\
	else if constexpr (std::same_as<ScalarType, f64>)											\
	{																							\
		status = operation(d)(__VA_ARGS__);														\
	}																							\
	else if constexpr (std::same_as<ScalarType, c32>)											\
	{																							\
		status = operation(c)(__VA_ARGS__);														\
	}																							\
	else if constexpr (std::same_as<ScalarType, c64>)											\
	{																							\
		status = operation(z)(__VA_ARGS__);														\
	}																							\
	AssertE(status == SPARSE_STATUS_SUCCESS, MessageTag::Math,									\
			Format("{}: {}.", name, MKL::SparseStatusMessages[status]));	

#define MKL_VML_OPERATION(operation, ...)														\
	if constexpr (std::same_as<ScalarType, f32>)												\
	{																							\
		operation(s)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, f64>)											\
	{																							\
		operation(d)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, c32>)											\
	{																							\
		operation(c)(__VA_ARGS__);																\
	}																							\
	else if constexpr (std::same_as<ScalarType, c64>)											\
	{																							\
		operation(z)(__VA_ARGS__);																\
	}