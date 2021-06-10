#pragma once

#include "Utils/Aliases.h"

namespace CESDSOL
{
	template<typename MatrixType>
	class Preconditioner
	{
	public:
		struct OutputInfo
		{
			bool success = true;
		};

		virtual uptr<OutputInfo> Apply(MatrixType& matrix) const noexcept = 0;
		virtual ~Preconditioner() = default;
	};
}