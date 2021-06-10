#pragma once

#include <complex>
#include <cstdint>
#include <memory>
#include <optional>

namespace CESDSOL
{
	using i8 = std::int8_t;
	using i16 = std::int16_t;
	using i32 = std::int32_t;
	using i64 = std::int64_t;

	using u8 = std::uint8_t;
	using u16 = std::uint16_t;
	using u32 = std::uint32_t;
	using u64 = std::uint64_t;

	using f32 = float;
	using f64 = double;

	using c32 = std::complex<f32>;
	using c64 = std::complex<f64>;

	template<typename T>
	using uptr = std::unique_ptr<T>;

	template<typename T>
	using sptr = std::shared_ptr<T>;

	template<typename T>
	using ref = std::reference_wrapper<T>;

	template<typename T>
	using cref = std::reference_wrapper<const T>;

	template<typename T>
	using opt = std::optional<T>;

	template<typename T>
	using optref = opt<ref<T>>;

	template<typename T>
	using optcref = opt<cref<T>>;
}
