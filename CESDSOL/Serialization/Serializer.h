#pragma once

#include "Grid/Grid.h"
#include "Serialization/DataToLoad.h"
#include "Serialization/GridCache.h"
#include "Utils/Utils.h"

#include <filesystem>
#include <fstream>
#include <string>

namespace CESDSOL::Serializer
{
	namespace fs = std::filesystem;
	
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

	template<typename T>
	void Write(std::ofstream& stream, const T& value) noexcept
	{
		stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
	}

	template<typename T>
	T Read(std::ifstream& stream) noexcept
	{
		T result;
		stream.read(reinterpret_cast<char*>(&result), sizeof(T));
		return result;
	}

	template<typename T>
	void WriteData(std::ofstream& stream, const T& object) noexcept
	{
		stream.write(reinterpret_cast<const char*>(object.data()), object.size() * sizeof(T::value_type));
	}

	[[nodiscard]] Array<uint8_t> ReadData(std::ifstream& stream, size_t dataSize) noexcept
	{
		Array<uint8_t> result(dataSize);
		stream.read(reinterpret_cast<char*>(result.data()), dataSize);
		return result;
	}

	constexpr uint32_t CurrentVersion = 1;
	constexpr uint32_t MagicNumber = 0x44534543;
	
	template<typename ProblemType>
	void Save(ProblemType& problem, const std::string& path) noexcept
	{
		auto stream = OpenFileForWrite(path);
		Write(stream, MagicNumber);
		Write(stream, CurrentVersion);
		const auto gridData = problem.GetGrid().Save();
		const auto problemData = problem.Save();
		Write(stream, gridData.type);
		Write(stream, gridData.header.size()); // Grid header size.
		Write(stream, gridData.data.size()); // Grid data size.
		WriteData(stream, gridData.header);
		WriteData(stream, gridData.data);
		Write(stream, problemData.type);
		Write(stream, sizeof(problemData.header)); // Problem header size.
		Write(stream, problemData.data.size()); // Problem data size.
		Write(stream, problemData.header);
		WriteData(stream, problemData.data);
	}	

	template<typename ProblemType>
	void Load(ProblemType& problem, const std::string& path, DataToLoad dataToLoad = DataToLoad::Everything) noexcept
	{
		auto stream = OpenFileForRead(path);
		AssertE(Read<uint32_t>(stream) == MagicNumber, MessageTag::Serialization, 
			Format("File {} has unknown format!", path));
		AssertE(Read<uint32_t>(stream) <= CurrentVersion, MessageTag::Serialization, 
			Format("File {} has unsupported version!", path));
		const auto gridType = Read<GridTypeData>(stream);
		AssertE(gridType.dimension == ProblemType::Dimension, MessageTag::Serialization, 
			"Grid dimension is inconsistent with problem dimension!");
		AssertE(static_cast<DataType>(gridType.dataType) == Serializer::ToDataType<ProblemType::CoordinateType>(), MessageTag::Serialization,
			"Grid data type is inconsistent with problem coordinate type!");
		const auto gridHeaderSize = Read<size_t>(stream);
		const auto gridDataSize = Read<size_t>(stream);
		const auto gridHeader = ReadData(stream, gridHeaderSize);
		const auto gridData = ReadData(stream, gridDataSize);
		const auto grid = GridCache::Instance().Load(gridType, gridHeader, gridData);
		AssertE(grid != nullptr, MessageTag::Serialization, "Cannot load grid!");
		AssertE(Read<uint32_t>(stream) == static_cast<uint32_t>(ProblemType::SerializerProblemType), MessageTag::Serialization, 
			"Trying to read problem data of inconsistent type!");
		const auto problemHeaderSize = Read<size_t>(stream);
		const auto problemDataSize = Read<size_t>(stream);
		const auto header = ReadData(stream, problemHeaderSize);
		const auto data = ReadData(stream, problemDataSize);
		problem.Load(header, data, std::static_pointer_cast<typename ProblemType::GridType>(grid), dataToLoad);
	}
}
