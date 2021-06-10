#pragma once

#include <cstdint>

namespace CESDSOL::Serializer
{
	struct GridTypeData
	{
		uint64_t type;
		uint64_t dimension;
		uint64_t dataType;

		auto operator<=>(const GridTypeData& other) const noexcept = default;
	};
}

namespace std
{
	template<>
	class hash<CESDSOL::Serializer::GridTypeData>
	{
	public:
		[[nodiscard]] size_t operator()(const CESDSOL::Serializer::GridTypeData& key) const noexcept
		{
			const auto hasher = std::hash<uint64_t>();
			auto result = hasher(key.type);
			result = HashCombine(result, key.dimension, hasher);
			result = HashCombine(result, key.dataType, hasher);
			return result;
		}
	};
}