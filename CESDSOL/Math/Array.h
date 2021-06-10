#pragma once

#include "Utils/Utils.h"

#include <cstddef>
#include <memory>
#include <iostream>

namespace CESDSOL
{
	template<typename ScalarType>
	class Array
	{
	protected:
		ScalarType* values = nullptr;
		size_t count = 0;

		void Free() noexcept
		{
			delete values;
		}

	public:
		using value_type = ScalarType;
		using size_type = size_t;

		[[nodiscard]] constexpr size_t size() const noexcept
		{
			return count;
		}

		[[nodiscard]] constexpr ScalarType* data() const noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr ScalarType* data() noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr ScalarType* begin() noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr ScalarType* end() noexcept
		{
			return values + size();
		}

		[[nodiscard]] constexpr ScalarType* begin() const noexcept
		{
			return values;
		}

		[[nodiscard]] constexpr ScalarType* end() const noexcept
		{
			return values + size();
		}

		constexpr Array() noexcept = default;

		constexpr explicit Array(size_t aCount, const ScalarType* source = nullptr) noexcept :
			count(aCount)
		{
			if (aCount > 0)
			{
				values = static_cast<ScalarType*>(::operator new(aCount * sizeof(ScalarType), std::nothrow));
				AssertE(values != nullptr, MessageTag::Math, Format("Memory allocation for array of size {} failed.", count));

				if (source == nullptr)
				{
					std::uninitialized_value_construct_n(values, count);
				}
				else
				{
					std::uninitialized_copy(source, source + count, values);
				}
			}
		}

		constexpr Array(ScalarType* data, size_t aCount) noexcept :
			values(data),
			count(aCount)
		{}

		constexpr Array(std::initializer_list<ScalarType> data) noexcept :
			Array(data.size(), data.begin())
		{}

		constexpr Array(const ScalarType& value, size_t count) noexcept :
			Array(count)
		{
			std::fill_n(begin(), count, value);
		}

		constexpr Array(Array<ScalarType>&& other) noexcept:
			Array(other.values, other.size())
		{
			other.values = nullptr;
			other.count = 0;
		}

		constexpr Array(const Array<ScalarType>& other) noexcept :
			Array(other.size(), other.begin())
		{}

		constexpr Array& operator=(Array<ScalarType>&& other) noexcept
		{
			if (&other != this)
			{
				std::swap(values, other.values);
				std::swap(count, other.count);
			}
			return *this;
		}

		constexpr Array& operator=(const Array<ScalarType>& other) noexcept
		{
			if (&other != this)
			{
				if (other.size() != size())
				{
					Free();
					if (other.size() > 0)
					{
						values = static_cast<ScalarType*>(::operator new(other.size() * sizeof(ScalarType), std::nothrow));
						AssertE(values != nullptr, MessageTag::Math, Format("Memory allocation for array of size {} failed.", count));

						count = other.size();
						std::uninitialized_copy(other.begin(), other.end(), values);
					}
					else
					{
						values = nullptr;
					}
				}
				else
				{
					std::copy(other.begin(), other.end(), values);
				}
			}
			return *this;
		}

		~Array()
		{
			Free();
		}

		void Release() noexcept
		{
			values = nullptr;
		}

		[[nodiscard]] constexpr const ScalarType& operator[](size_t index) const noexcept
		{
			return values[index];
		}

		[[nodiscard]] constexpr ScalarType& operator[](size_t index) noexcept
		{
			return values[index];
		}

	public:
		[[nodiscard]] constexpr const ScalarType& At(size_t index) const noexcept
		{
			AssertE(index < count, MessageTag::Math, Format("Attempt to read element {} in array of size {}.", index, count));
			return values[index];
		}

		[[nodiscard]] constexpr ScalarType& At(size_t index) noexcept
		{
			AssertE(index < count, MessageTag::Math, Format("Attempt to write element {} to array of size {}.", index, count));
			return values[index];
		}

		[[nodiscard]] constexpr ScalarType Front() const noexcept
		{
			return values[0];
		}

		[[nodiscard]] constexpr ScalarType Back() const noexcept
		{
			return values[size() - 1];
		}
	};
}
