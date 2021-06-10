#pragma once

#include "Utils/Aliases.h"

#include <string>

namespace CESDSOL
{
	enum class MessageTag
	{
		Core,
		Math,
		Problem,
		Input,
		Output,
		Discretization,
		Jacobian,
		LinearSolver,
		NonlinearSolver,
		Preconditioner,
		ParametricSweeper,
		LineSearcher,
		Serialization,
		TransientSolver,
		Count
	};

	enum class MessageType
	{
		Info,
		Warning,
		Error,
		Debug,
		Count
	};

	enum class MessagePriority
	{
		Low,
		Medium,
		High,
		Critical,
		Count
	};

	class LogFormatter
	{
	public:
		[[nodiscard]] virtual std::string Format(MessageType type, MessagePriority priority, MessageTag tag, const std::string& message) const noexcept = 0;
		virtual ~LogFormatter() = default;
	};

	class LogWriter
	{
	public:
		[[nodiscard]] virtual bool CanWrite(MessageType type, MessagePriority priority, MessageTag tag) const noexcept = 0;
		virtual void Write(const std::string& message) const noexcept = 0;
		virtual ~LogWriter() = default;
	};

	class DefaultLogFormatter final : public LogFormatter
	{
	public:
		[[nodiscard]] std::string Format(MessageType type, MessagePriority priority, MessageTag tag, const std::string& message) const noexcept override;
	};

	class DefaultLogWriter final : public LogWriter
	{
	public:
		[[nodiscard]] bool CanWrite(MessageType type, MessagePriority priority, MessageTag tag) const noexcept override;
		void Write(const std::string& message) const noexcept override;
	};

	namespace Logger
	{		
		void Log(MessageType type, MessagePriority priority, MessageTag tag, const std::string& message) noexcept;
		void SetFormatter(uptr<LogFormatter> aFormatter) noexcept;
		void SetWriter(uptr<LogWriter> aWriter) noexcept;
	}
}
