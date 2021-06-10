#pragma once

#include "Grid/BaseGrid.h"
#include "Math/Array.h"
#include "Serialization/GridTypeData.h"
#include "Utils/Singleton.h"

namespace CESDSOL::Serializer
{
	class GridCache : public Singleton<GridCache>
	{
	public:
		void Register(const GridTypeData& typeData, const std::function<std::shared_ptr<BaseGrid>(const Array<uint8_t>&, const Array<uint8_t>&)>& loader) noexcept;

		[[nodiscard]] std::shared_ptr<BaseGrid> Load(const GridTypeData& type, const Array<uint8_t>& header, const Array<uint8_t>& data) noexcept;

	private:
		std::unordered_map<GridTypeData, std::function<std::shared_ptr<BaseGrid>(const Array<uint8_t>&, const Array<uint8_t>&)>> entries;
	};
}

#define RegisterGridLoader(typeData, loader) inline static int gridRegistrator = (CESDSOL::Serializer::GridCache::Instance().Register(typeData, loader), 0)