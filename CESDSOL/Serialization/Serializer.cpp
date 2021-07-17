#include "Serialization/Serializer.h"

namespace CESDSOL::Serializer
{
	[[nodiscard]] std::ofstream OpenFileForWrite(const fs::path& path) noexcept
	{
		const auto parentPath = path.parent_path();
		create_directories(parentPath);
		std::ofstream result(path, std::ios::binary);
		const bool isGood = result.good();
		return result;
	}

	[[nodiscard]] std::ifstream OpenFileForRead(const fs::path& path) noexcept
	{
		std::ifstream result(path, std::ios::binary);
		AssertE(result.good(), MessageTag::Serialization, Format("Failed to open file {} for read!", path.string()));
		return result;
	}

	[[nodiscard]] Array<uint8_t> ReadData(std::ifstream& stream, size_t dataSize) noexcept
	{
		Array<uint8_t> result(dataSize);
		stream.read(reinterpret_cast<char*>(result.data()), dataSize);
		return result;
	}
}