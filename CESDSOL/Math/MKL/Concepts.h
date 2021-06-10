#pragma once

#include "Utils/Aliases.h"

#include <concepts>

namespace CESDSOL::MKL
{
	template<typename T>
	concept MKLRealScalar = std::same_as<T, f32> || std::same_as<T, f64>;

	template<typename T>
	concept MKLScalar = MKLRealScalar<T> || std::same_as<T, c32> || std::same_as<T, c64>;
}