#include "Math/MKL/ILUT.h"

#include "ILU0.h"

namespace CESDSOL::MKL
{
	uptr<ILUT::OutputInfo> ILUT::Apply(CSRMatrix<double>& matrix) const noexcept
	{
		const int size = matrix.RowCount();
		int error;
		const size_t newNonzeroCount = (2 * maxfil + 1) * size - maxfil * (maxfil + 1) + 1;
		CSRMatrix<double> preconditioned(size, size, newNonzeroCount);
		
		dcsrilut(&size, matrix.GetValues().data(), matrix.GetRowCounts().data(), matrix.GetColumnIndices().data(),
			preconditioned.GetValues().data(), preconditioned.GetRowCounts().data(), preconditioned.GetColumnIndices().data(), 
			&tolerance, &maxfil, intParameters.data(), doubleParameters.data(), &error);

		switch (error)
		{
		case 0:
		case 101:
		case 102:
		case 103:
		case 104:
			matrix = preconditioned;
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"ILUT preconditioner was successfully calculated.");
			return std::make_unique<OutputInfo>(true);

		case -101:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILUT preconditioner calculation: the number of elements in some matrix row is equal to or less than 0.");
			break;

		case -102:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILUT preconditioner calculation: computed diagonal element is less than the product of the given tolerance and the current matrix row norm.");
			break;

		case -103:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILUT preconditioner calculation: incorrect CSR matrix element ordering.");
			break;

		case -104:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILUT preconditioner calculation: memory is insufficient for the internal work array.");
			break;

		case -105:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILUT preconditioner calculation: negative value of maximum fill-in parameter.");
			break;

		case -107:
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::Preconditioner,
				"Error in ILU0 preconditioner calculation: CSR matrix contains invalid column index.");
			break;
		}
		return std::make_unique<OutputInfo>(false);
	}

	ILUT::ILUT() noexcept
		: intParameters()
		, doubleParameters()
	{
		intParameters[5] = 0;
		intParameters[6] = 0;
		intParameters[30] = 1;
		doubleParameters[30] = 1e-10;
	}

	void ILUT::SetNormalizeZeroDiagonal(bool value) noexcept
	{
		intParameters[30] = static_cast<int>(value);
	}

	void ILUT::SetZeroDiagonalNormalizer(double value) noexcept
	{
		doubleParameters[30] = value;
	}
}
