#pragma once

#include "Math/MultiLevelArray.h"

#include <map>

namespace CESDSOL
{
	template<typename TimeType = double, typename FieldType = double>
	class TimeCache
	{
	public:
		void Add(TimeType time, const TwoLevelArray<FieldType>& solution) noexcept
		{
			if (!cache.empty() && cache.rbegin()->first > time)
			{
				ClearSince(time);
			}
			cache.insert_or_assign(time, solution);
		}

		void Clear() noexcept
		{
			cache.clear();
		}

		void ClearSince(TimeType time) noexcept
		{
			auto it = cache.lower_bound(time);
			cache.erase(it);
		}

		[[nodiscard]] size_t CachedCount() const noexcept
		{
			return cache.size();
		}

		[[nodiscard]] const std::map<TimeType, TwoLevelArray<FieldType>>& GetData() const noexcept
		{
			return cache;
		}

		[[nodiscard]] optcref<TwoLevelArray<FieldType>> Get(TimeType time) const noexcept
		{
			if (auto it = cache.find(time); it != cache.end())
			{
				return std::cref(it->second);
			}
			return std::nullopt;
		}

		[[nodiscard]] optcref<TwoLevelArray<FieldType>> GetNext(TimeType time) const noexcept
		{
			if (auto it = cache.lower_bound(time); it != cache.end())
			{
				return std::cref(it->second);
			}
			return std::nullopt;
		}

		[[nodiscard]] optcref<TwoLevelArray<FieldType>> GetPrevious(TimeType time) const noexcept
		{
			auto it = cache.lower_bound(time);
			if (it != cache.begin())
			{
				if (it == cache.end() || it->first != time)
				{
					return std::cref((--it)->second);
				}
				else
				{
					return std::cref(it->second);
				}
			}			
			return std::nullopt;
		}

	private:		
		std::map<TimeType, TwoLevelArray<FieldType>> cache;
	};
}