#pragma once
namespace CESDSOL
{
	template<typename IndexType, typename BodyType>
	void ParallelFor(IndexType startIndex, IndexType endIndex, BodyType&& body) noexcept
	{
		for (IndexType index = startIndex; index < endIndex; ++index)
		{
			body();
		}
	}

	template<typename BodyType>
	void ParallelBlock(BodyType&& body) noexcept
	{
		body();
	}

	template<typename IndexType, typename BodyType>
	void ForInParallelBlock(IndexType startIndex, IndexType endIndex, BodyType&& body) noexcept
	{
		for (IndexType index = startIndex; index < endIndex; ++index)
		{
			body();
		}
	}
}