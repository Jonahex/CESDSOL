#pragma once

#include "Serialization/GridCache.h"

namespace CESDSOL::Serializer
{
	void GridCache::Register(const GridTypeData& typeData, const std::function<std::shared_ptr<BaseGrid>(const Array<uint8_t>&, const Array<uint8_t>&)>& loader) noexcept
	{
		entries.insert_or_assign(typeData, loader);
	}

	[[nodiscard]] std::shared_ptr<BaseGrid> GridCache::Load(const GridTypeData& type, const Array<uint8_t>& header, const Array<uint8_t>& data) noexcept
	{
		const auto it = entries.find(type);
		if (it != entries.end())
		{
			return it->second(header, data);
		}
		return nullptr;
	}
}