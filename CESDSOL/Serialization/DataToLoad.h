#pragma once

#include "Utils/Utils.h"

namespace CESDSOL::Serializer
{
	enum class DataToLoad : uint64_t
	{
		Variables = 1 << 0,
		Parameters = 1 << 1,
		Time = 1 << 2,
		Everything = 7,
	};
	MakeFlag(DataToLoad)
}
