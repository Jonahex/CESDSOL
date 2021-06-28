#pragma once

namespace CESDSOL
{
	template<typename BodyType>
	void ParallelFor(int64_t startIndex, int64_t endIndex, BodyType&& body) noexcept
	{
#pragma omp parallel for
		for (int64_t index = startIndex; index < endIndex; ++index)
		{
			body(index);
		}
	}

	template<typename BodyType>
	void ParallelBlock(BodyType&& body) noexcept
	{
#pragma omp parallel
		body();
	}

	template<typename BodyType>
	void ForInParallelBlock(int64_t startIndex, int64_t endIndex, BodyType&& body) noexcept
	{
#pragma omp for
		for (int64_t index = startIndex; index < endIndex; ++index)
		{
			body(index);
		}
	}
}