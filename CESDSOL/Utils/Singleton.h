#pragma once

namespace CESDSOL
{
	template<typename T>
	class Singleton
	{
	public:
		[[nodiscard]] static T& Instance() noexcept
		{
			static T instance;
			return instance;
		}
	};
}