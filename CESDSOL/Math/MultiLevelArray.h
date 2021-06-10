#pragma once

#include "Math/Array.h"

#include <algorithm>
#include <numeric>
#include <span>

namespace CESDSOL
{
	template<size_t LevelCount> requires (LevelCount > 1)
	class LevelStructure final
		: public Array<std::pair<size_t, LevelStructure<LevelCount - 1>>>
	{
	public:
		constexpr LevelStructure() noexcept = default;

		LevelStructure(size_t first, const LevelStructure<LevelCount - 1>& second) :
			Array<std::pair<size_t, LevelStructure<LevelCount - 1>>>({ {first, second} })
		{}

		LevelStructure(std::initializer_list<std::pair<size_t, LevelStructure<LevelCount - 1>>> elements) :
			Array<std::pair<size_t, LevelStructure<LevelCount - 1>>>(elements)
		{}

		LevelStructure(const Array<std::pair<size_t, LevelStructure<LevelCount - 1>>>& elements) :
			Array<std::pair<size_t, LevelStructure<LevelCount - 1>>>(elements)
		{}

		[[nodiscard]] size_t ElementCount() const noexcept
		{
			return std::accumulate(this->begin(), this->end(), 0, 
				[](auto sum, const auto& element) {return sum + element.first * element.second.ElementCount(); });
		}
	};

	template<>
	class LevelStructure<2> final 
		: public Array<std::pair<size_t, size_t>>
	{
	public:
		constexpr LevelStructure() noexcept = default;

		LevelStructure(size_t first, size_t second) :
			Array<std::pair<size_t, size_t>>({ {first, second} })
		{}

		LevelStructure(std::initializer_list<std::pair<size_t, size_t>> elements) :
			Array<std::pair<size_t, size_t>>(elements)
		{}

		LevelStructure(const Array<std::pair<size_t, size_t>> & elements) :
			Array<std::pair<size_t, size_t>>(elements)
		{}

		[[nodiscard]] size_t ElementCount() const noexcept
		{
			return std::accumulate(this->begin(), this->end(), 0,
				[](auto sum, const auto& element) {return sum + element.first * element.second; });
		}
	};

	template<size_t LevelCount>
	[[nodiscard]] size_t SubspanCount(const LevelStructure<LevelCount>& structure) noexcept
	{
		return std::accumulate(structure.begin(), structure.end(), 0,
			[](auto sum, const auto& element) {return sum + element.first; });
	}

	template<typename ScalarType, size_t LevelCount>
	class MultiLevelSpan 
		: public Array<MultiLevelSpan<ScalarType, LevelCount - 1>>
	{
	public:
		constexpr MultiLevelSpan() noexcept = default;

		constexpr MultiLevelSpan(ScalarType* data, const LevelStructure<LevelCount>& structure):
			Array<MultiLevelSpan<ScalarType, LevelCount - 1>>(SubspanCount(structure))
		{
			size_t index = 0;
			for (const auto& element : structure)
			{
				auto shift = element.second.ElementCount();
				for (size_t i = 0; i < element.first; i++)
				{
					this->operator[](index++) = MultiLevelSpan<ScalarType, LevelCount - 1>(data, element.second);
					data += shift;
				}
			}
		}
	}; 

	template<typename ScalarType>
	class MultiLevelSpan<ScalarType, 2> 
		: public Array<std::span<ScalarType>>
	{
	public:
		constexpr MultiLevelSpan() noexcept = default;

		constexpr MultiLevelSpan(ScalarType* data, const LevelStructure<2>& structure) :
			Array<std::span<ScalarType>>(SubspanCount(structure))
		{
			size_t index = 0;
			for (const auto& element : structure)
			{
				for (size_t i = 0; i < element.first; i++)
				{
					this->operator[](index++) = std::span(data, element.second);
					data += element.second;
				}
			}
		}
	};

	template<typename ScalarType, size_t LevelCount> requires (LevelCount > 1)
	class MultiLevelArray
		: private Array<ScalarType>
		, public MultiLevelSpan<ScalarType, LevelCount>
	{
	public:
		using MultiLevelSpan<ScalarType, LevelCount>::value_type;
		using MultiLevelSpan<ScalarType, LevelCount>::size_type;
		using MultiLevelSpan<ScalarType, LevelCount>::begin;
		using MultiLevelSpan<ScalarType, LevelCount>::end;
		using MultiLevelSpan<ScalarType, LevelCount>::size;
		using MultiLevelSpan<ScalarType, LevelCount>::data;
		using MultiLevelSpan<ScalarType, LevelCount>::operator[];

		constexpr MultiLevelArray() noexcept = default;

		constexpr MultiLevelArray(const LevelStructure<LevelCount>& structure) noexcept:
			Array<ScalarType>(structure.ElementCount()),
			MultiLevelSpan<ScalarType, LevelCount>(Array<ScalarType>::data(), structure)
		{}

		[[nodiscard]] constexpr const Array<ScalarType>& Flatten() const noexcept
		{
			return *this;
		}

		constexpr Array<ScalarType>& Flatten() noexcept
		{
			return *this;
		}

	private:
		void ShiftLevels(std::intptr_t shift) noexcept
		{
			if constexpr (LevelCount == 2)
			{
				for (auto& element : *this)
				{
					element = std::span<ScalarType>(element.data() + shift, element.size());
				}
			}
			else
			{
				for (auto& element : *this)
				{
					element.ShiftLevels(shift);
				}
			}
		}

	public:
		constexpr MultiLevelArray<ScalarType, LevelCount>& operator=(const Array<ScalarType>& other) noexcept
		{
			Flatten() = other;
			return *this;
		}

		constexpr MultiLevelArray<ScalarType, LevelCount>& operator=(Array<ScalarType>&& other) noexcept
		{
			ShiftLevels(other.data() - Flatten().data());
			Flatten() = std::move(other);
			return *this;
		}
	};

	template<typename ScalarType>
	using TwoLevelArray = MultiLevelArray<ScalarType, 2>;

	template<typename ScalarType>
	using ThreeLevelArray = MultiLevelArray<ScalarType, 3>;

	template<typename ScalarType>
	using FourLevelArray = MultiLevelArray<ScalarType, 4>;

	template<typename ScalarType, size_t LevelCount>
	[[nodiscard]] LevelStructure<LevelCount> ApproximateStructure(const MultiLevelSpan<ScalarType, LevelCount>& source) noexcept
	{
		if constexpr (LevelCount == 2)
		{
			Array<std::pair<size_t, size_t>> structure{ source.size() };
			std::transform(source.begin(), source.end(), structure.begin(), [](const auto& x) -> std::pair<size_t, size_t> {return { 1, x.size() }; });
			return structure;
		}
		else
		{
			Array<std::pair<size_t, LevelStructure<LevelCount - 1>>> structure{ source.size() };
			std::transform(source.begin(), source.end(), structure.begin(), [](const auto& x) -> std::pair<size_t, LevelStructure<LevelCount - 1>> {return { 1, ApproximateStructure(x) }; });
			return structure;
		}
	}

	template<typename ScalarType, size_t LevelCount, typename TransformType>
	[[nodiscard]] MultiLevelArray<ScalarType, LevelCount> MimicStructure(const MultiLevelArray<ScalarType, LevelCount>& source) noexcept
	{
		return MultiLevelArray<ScalarType, LevelCount>(ApproximateStructure(source));
	}
}