#pragma once

#include <concepts>
#include <cstdint>

namespace CESDSOL::Serializer
{
	enum class DataType
	{
		Int8,
		Int16,
		Int32,
		Int64,
		UInt8,
		UInt16,
		UInt32,
		UInt64,
		Float16,
		Float32,
		Float64,
		Undefined
	};

	template<typename T>
	[[nodiscard]] constexpr DataType ToDataType() noexcept
	{
		if constexpr (std::same_as<T, std::int8_t>) return DataType::Int8;
		else if constexpr (std::same_as<T, std::int16_t>) return DataType::Int16;
		else if constexpr (std::same_as<T, std::int32_t>) return DataType::Int32;
		else if constexpr (std::same_as<T, std::int64_t>) return DataType::Int64;
		else if constexpr (std::same_as<T, std::uint8_t>) return DataType::UInt8;
		else if constexpr (std::same_as<T, std::uint16_t>) return DataType::UInt16;
		else if constexpr (std::same_as<T, std::uint32_t>) return DataType::UInt32;
		else if constexpr (std::same_as<T, std::uint64_t>) return DataType::UInt64;
		else if constexpr (std::same_as<T, float>) return DataType::Float32;
		else if constexpr (std::same_as<T, double>) return DataType::Float64;
		return DataType::Undefined;
	}
}