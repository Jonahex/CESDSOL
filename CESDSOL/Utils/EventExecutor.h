#pragma once

namespace CESDSOL
{
	template<typename ExecutorType, typename ActionTargetType, typename Event>
	class EventExecutor
	{
	public:
		void AddAction(Event event, const std::function<void(ActionTargetType&)>& action) noexcept
		{
			events[static_cast<size_t>(event)].push_back(action);
		}

		void RemoveAction(Event event, const std::function<void(ActionTargetType&)>& action) noexcept
		{
			auto& actions = events[static_cast<size_t>(event)];
			actions.erase(std::remove(actions.begin(), actions.end(), action));
		}

	protected:
		void ApplyActions(Event event, ActionTargetType& target) const noexcept
		{
			for (const auto& action : events[static_cast<size_t>(event)])
			{
				action(target);
			}
		}

	private:
		std::array<std::vector<std::function<void(ActionTargetType&)>>, static_cast<size_t>(Event::Count)> events;
	};
}