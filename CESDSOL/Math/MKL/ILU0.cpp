#include "Math/MKL/ILU0.h"

namespace CESDSOL::MKL
{
	uptr<ILU0::OutputInfo> ILU0::Apply(CSRMatrix<double>& matrix) const noexcept
	{
		const int size = matrix.RowCount();
		int error;
		Array<double> preconditioner(matrix.NonZeroCount());
		dcsrilu0(&size, matrix.GetValues().data(), matrix.GetRowCounts().data(), matrix.GetColumnIndices().data(),
			preconditioner.data(), intParameters.data(), doubleParameters.data(), &error);

		switch (error)
		{
			case 0:
				matrix.ReplaceValues(preconditioner);
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner, 
					"ILU0 preconditioner was successfully calculated.");
				return std::make_unique<OutputInfo>(true);

			case -101:
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
					"Error in ILU0 preconditioner calculation: at least one diagonal element is omitted from the matrix.");
				break;

			case -102:
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
					"Error in ILU0 preconditioner calculation: the matrix contains a diagonal element with the value of zero.");
				break;

			case -103:
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
					"Error in ILU0 preconditioner calculation: the matrix contains a diagonal element which is too small.");
				break;

			case -104:
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
					"Error in ILU0 preconditioner calculation: memory is insufficient for the internal work array.");
				break;

			case -106:
				Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
					"Error in ILU0 preconditioner calculation: matrix column indices are not in the ascending order.");
				break;
		}
		return std::make_unique<OutputInfo>(false);
	}

	ILU0::ILU0() noexcept
		: intParameters()
		, doubleParameters()
	{
		intParameters[5] = 0;
		intParameters[30] = 1;
	}
	
	void ILU0::SetNormalizeZeroDiagonal(bool value) noexcept
	{
		intParameters[30] = static_cast<int>(value);
	}
	
	void ILU0::SetZeroDiagonalThreshold(double value) noexcept
	{
		doubleParameters[30] = value;
	}
	
	void ILU0::SetZeroDiagonalNormalizer(double value) noexcept
	{
		doubleParameters[31] = value;
	}
}