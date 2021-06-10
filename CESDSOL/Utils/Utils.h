#pragma once

#include "Utils/Logger.h"
#include "Configuration.h"

#include <array>
#include <cstddef>
#include <format>
#include <functional>
#include <memory>
#include <type_traits>
#include <chrono>

#define _USE_MATH_DEFINES
#include <math.h>

#define MakeFlag(type) \
	type operator|(type left, type right) \
	{ \
		return static_cast<type>(static_cast<size_t>(left) | static_cast<size_t>(right)); \
	} \
	bool operator&(type left, type right) \
	{ \
		return static_cast<bool>(static_cast<size_t>(left) & static_cast<size_t>(right)); \
	}

#define MakeProperty(name, propertyName, type, defaultValue) \
	private: \
		type name = defaultValue; \
	public: \
		void Set##propertyName(type value) \
		{ \
			name = value; \
		} \
		[[nodiscard]] type Get##propertyName() const\
		{ \
			return name; \
		}

namespace CESDSOL
{
	template<typename Type> requires (std::is_floating_point_v<Type>)
	struct NiceFloat
	{
		Type value;
	};

	template<typename DurationType>
	concept Duration = requires (DurationType value)
	{
		std::chrono::duration(value);
	};

	template<Duration DurationType>
	struct NiceTime
	{
		std::chrono::microseconds::rep value;

		NiceTime(const DurationType& duration)
			: value(std::chrono::duration_cast<std::chrono::microseconds>(duration).count())
		{}
	};

	template<typename Type>
	[[nodiscard]] auto Nicify(const Type& value) noexcept
	{
		if constexpr (std::is_floating_point_v<Type>)
		{
			return NiceFloat{ value };
		}
		else if constexpr (Duration<Type>)
		{
			return NiceTime(value);
		}
		else return value;
	}

	template<typename... Args>
	[[nodiscard]] std::string Format(std::string_view format, const Args&... args) noexcept
	{
		return std::format(format, Nicify(args)...);
	}

	template<typename OutputIt, typename... Args>
	[[nodiscard]] OutputIt FormatTo(OutputIt out, std::string_view format, const Args&... args) noexcept
	{
		return std::format_to(out, format, Nicify(args)...);
	}
}

template<typename FloatType, class CharType>
struct std::formatter<CESDSOL::NiceFloat<FloatType>, CharType>
{
	auto parse(std::format_parse_context& ctx)
	{
		return ctx.end();
	}

	template<class FormatContext>
	auto format(CESDSOL::NiceFloat<FloatType> fp, FormatContext& ctx)
	{
		std::string formatted;
		if (fp.value == 0)
		{
			formatted = "0";
		}
		const auto power = static_cast<int>(std::floor(std::log10(std::abs(fp.value))));
		if (power < 4 && power > -4)
		{
			formatted = std::format("{:.4f}", fp.value);
		}
		else
		{
			formatted = std::format("{:.4f}e{}", fp.value / std::pow(10, power), power);
		}
		ctx.advance_to(std::copy(formatted.begin(), formatted.end(), ctx.out()));
		return ctx.out();
	}
};


template<typename DurationType, class CharType>
struct std::formatter<CESDSOL::NiceTime<DurationType>, CharType>
{
	auto parse(std::format_parse_context& ctx)
	{
		return ctx.end();
	}

	template<class FormatContext>
	auto format(CESDSOL::NiceTime<DurationType> time, FormatContext& ctx)
	{
		std::string formatted;
		if (time.value < 1000)
		{
			formatted = CESDSOL::Format("{} mcs", time.value);
		}
		else if (time.value < 1000000)
		{
			formatted = CESDSOL::Format("{} ms", time.value / 1000.);
		}
		else if (time.value < 60000000)
		{
			formatted = CESDSOL::Format("{} s", time.value / 1000000.);
		}
		else if (time.value < 3600000000)
		{
			formatted = CESDSOL::Format("{} m {} s", time.value / 60000000, (time.value % 60000000) / 1000000.);
		}
		else if (time.value < 86400000000)
		{
			formatted = CESDSOL::Format("{} h {} m {} s", time.value / 3600000000, (time.value % 3600000000) / 60000000, (time.value % 60000000) / 1000000.);
		}
		else
		{
			formatted = CESDSOL::Format("{} d {} h {} m {} s", time.value / 86400000000, (time.value % 86400000000) / 3600000000, (time.value % 3600000000) / 60000000, (time.value % 60000000) / 1000000.);
		}
		ctx.advance_to(std::copy(formatted.begin(), formatted.end(), ctx.out()));
		return ctx.out();
	}
};

namespace CESDSOL
{
	const double GoldenRatio = (1 + std::sqrt(5)) / 2;
	constexpr double Pi = M_PI;
	constexpr double E = M_E;

	template<typename T>
	[[nodiscard]] T csc(T x) noexcept
	{
		return 1 / sin(x);
	}

	template<typename T>
	[[nodiscard]] T sec(T x) noexcept
	{
		return 1 / cos(x);
	}

	template<typename T>
	[[nodiscard]] T cot(T x) noexcept
	{
		return 1 / tan(x);
	}

	template<typename... Types>
	[[nodiscard]] std::array<std::common_type_t<Types...>, sizeof...(Types)> MakeArray(Types&&... elements) noexcept
	{
		return { std::forward<Types>(elements)... };
	}

#define AssertE(condition, tag, message) \
	if (!(condition)) \
	{ \
		Logger::Log(MessageType::Error, MessagePriority::Critical, tag, message); \
		std::terminate(); \
	}

#ifdef DebugMode
#define AssertD(condition, tag, message) AssertE(condition, tag, message)
#else
#define AssertD(condition, debug, message)
#endif

	template <typename T>
	[[nodiscard]] constexpr T Power(T value, size_t power) noexcept
	{
		T result = 1;
		for (size_t i = 0; i < power; i++)
		{
			result *= power;
		}
		return result;
	}

	template<typename Collection>
	[[nodiscard]] size_t CountNonZero(const Collection& collection) noexcept
	{
		size_t count = 0;
		for (const auto value : collection)
		{
			if (std::abs(value) > std::numeric_limits<typename Collection::value_type>::epsilon())
			{
				count++;
			}
		}
		return count;
	}
}

template<typename ValueType, typename HasherType>
[[nodiscard]] size_t HashCombine(size_t seed, const ValueType& value, const HasherType& hasher = std::hash<ValueType>()) noexcept
{
    return seed ^ (hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}

namespace std
{
	template<typename ValueType, size_t Dimension>
	class hash<std::array<ValueType, Dimension>>
	{
	public:
		[[nodiscard]] size_t operator()(const std::array<ValueType, Dimension>& key) const noexcept
		{
			const auto hasher = std::hash<ValueType>();
			auto result = hasher(key[0]);
			for (size_t i = 1; i < Dimension; i++)
			{
				result = HashCombine(result, key[i], hasher);
			}
			return result;
		}
	};
}