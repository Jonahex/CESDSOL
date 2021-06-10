#include "Utils/Logger.h"

#include "Utils/Utils.h"

#include <iostream>
#include <mutex>

namespace CESDSOL
{
	namespace
	{
		std::mutex logMutex;
		std::unique_ptr<LogFormatter> formatter = std::make_unique<DefaultLogFormatter>();
		std::unique_ptr<LogWriter> writer = std::make_unique<DefaultLogWriter>();

		constexpr std::array<const char*, static_cast<size_t>(MessageTag::Count)> MessageTagNames
		{
			"Core",
			"Math",
			"Problem",
			"Input",
			"Output",
			"Discretization",
			"Jacobian",
			"LinearSolver",
			"NonlinearSolver",
			"Preconditioner",
			"Traverser",
			"LineSearcher",
			"Serialization"
		};
		constexpr std::array<const char*, static_cast<size_t>(MessageType::Count)> MessageTypeNames
		{
		"Info",
		"Warning",
		"Error",
		"Debug"
		};
		constexpr std::array<const char*, static_cast<size_t>(MessagePriority::Count)> MessagePriorityNames
		{
			"Low",
			"Medium",
			"High",
			"Critical"
		};
	}

	void Logger::Log(MessageType type, MessagePriority priority, MessageTag tag, const std::string& message) noexcept
	{
		if (writer->CanWrite(type, priority, tag))
		{
			const auto guard = std::lock_guard(logMutex);
			writer->Write(formatter->Format(type, priority, tag, message));
		}
	}
	
	void Logger::SetFormatter(uptr<LogFormatter> aFormatter) noexcept
	{
		AssertE(aFormatter != nullptr, MessageTag::Core, "Trying to set invalid log formatter!");
		const auto guard = std::lock_guard(logMutex);
		formatter = std::move(aFormatter);
	}
	
	void Logger::SetWriter(uptr<LogWriter> aWriter) noexcept
	{
		AssertE(aWriter != nullptr, MessageTag::Core, "Trying to set invalid log writer!");
		const auto guard = std::lock_guard(logMutex);
		writer = std::move(aWriter);
	}

	[[nodiscard]] std::string DefaultLogFormatter::Format(MessageType type, MessagePriority priority, MessageTag tag, const std::string& message) const noexcept
	{
		if (type == MessageType::Info)
		{
			return message;
		}
		else
		{
			return CESDSOL::Format("[{}] [{}] {}", MessageTypeNames[static_cast<size_t>(type)], MessageTagNames[static_cast<size_t>(tag)], message.c_str());
		}
	}

	[[nodiscard]] bool DefaultLogWriter::CanWrite(MessageType type, MessagePriority priority, MessageTag tag) const noexcept
	{
#ifndef DebugMode
		if (type == MessageType::Debug)
		{
			return false;
		}
#endif
		MessagePriority priorityLimit = MessagePriority::Medium;
		if (tag == MessageTag::LineSearcher || tag == MessageTag::TransientSolver)
		{
			priorityLimit = MessagePriority::High;
		}
		return static_cast<size_t>(priority) >= static_cast<size_t>(priorityLimit);
	}

	void DefaultLogWriter::Write(const std::string& message) const noexcept
	{
		std::cout << message << '\n';
	}
}
