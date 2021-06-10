#pragma once

#include "Problem/BaseProblem.h"
#include "Problem/ExplicitTransientProblemDescriptor.h"
#include "Problem/TimeCache.h"

namespace CESDSOL
{
	template<size_t DimensionArg, template<typename> typename MatrixTypeArg = CSRMatrix, typename CoordinateTypeArg = double, typename FieldTypeArg = double>
	class ExplicitTransientProblem
		: public BaseProblem<ExplicitTransientProblemDescriptor<DimensionArg, CoordinateTypeArg, FieldTypeArg>, MatrixTypeArg>
	{
	private:
		using Self = ExplicitTransientProblem<DimensionArg, MatrixTypeArg, CoordinateTypeArg, FieldTypeArg>;
		using BaseType = BaseProblem<ExplicitTransientProblemDescriptor<DimensionArg, CoordinateTypeArg, FieldTypeArg>, MatrixTypeArg>;

	public:
		using typename BaseType::DescriptorType;
		using typename BaseType::CoordinateType;
		using typename BaseType::FieldType;
		using typename BaseType::GridType;
		using typename BaseType::DiscretizationType;
		using typename BaseType::CurrentGlobalValuesForVIEs;

		static constexpr size_t Dimension = BaseType::Dimension;

	private:
		using BaseType::isActualOnParameters;
		using BaseType::variables;

		using BaseType::BaseType;
		
		friend DescriptorType;

		CoordinateType time = 0;
		TimeCache<CoordinateType, FieldType> timeCache;
		
		[[nodiscard]] CurrentGlobalValuesForVIEs ConstructGlobalValuesForVIEs() const noexcept override
		{
			auto result = BaseType::ConstructGlobalValuesForVIEs();
			result.Time = time;
			return result;
		}

	public:
		[[nodiscard]] CoordinateType GetTime() const noexcept
		{
			return time;
		}

		void SetTime(CoordinateType aTime) noexcept
		{
			time = aTime;
			isActualOnParameters = false;
		}

		[[nodiscard]] const TimeCache<CoordinateType, FieldType>& GetCache() const noexcept
		{
			return timeCache;
		}

		[[nodiscard]] TimeCache<CoordinateType, FieldType>& GetCache() noexcept
		{
			return timeCache;
		}

		void CacheCurrent() noexcept
		{
			timeCache.Add(time, variables);
		}

		static constexpr Serializer::ProblemType SerializerProblemType = Serializer::ProblemType::TransientProblem;

		struct SerializedData
		{
			uint32_t type = static_cast<uint32_t>(SerializerProblemType);
			struct Header : BaseType::SerializedData::Header
			{
				uint64_t cachedCount;

				Header() noexcept = default;

				Header(const Self& problem) noexcept
					: BaseType::SerializedData::Header(problem)
					, cachedCount(problem.GetCache().CachedCount())
				{}
			} header;
			Array<uint8_t> data;
		};

	public:
		[[nodiscard]] SerializedData Save() noexcept
		{
			SerializedData result;
			result.header = SerializedData::Header(*this);

			const auto& cache = GetCache().GetData();

			const auto dataSize = sizeof(FieldType) *
				(this->GetDescriptor().ParameterCount() + (this->GetGrid().GetSize() * this->GetSerializedDataFieldCount() + this->GetSerializedDataVariableCount()) * cache.size())
				+ this->GetSerializedDataStringLength()
				+ sizeof(CoordinateType) * cache.size();

			result.data = Array<uint8_t>(dataSize);
			auto data = result.data.data();
			data = this->WriteSerializedDataParameters(data);
			for (const auto& [time, solution] : cache)
			{
				*reinterpret_cast<CoordinateType*>(data) = time;
				data += sizeof(CoordinateType);
				this->SetVariables(solution.Flatten());
				data = this->WriteSerializedDataFields(data);
				data = this->WriteSerializedDataLocalOutput(data);
				data = this->WriteSerializedDataGlobalOutput(data);
			}
			this->WriteSerializedDataStrings(data);

			return result;
		}

		void Load(const Array<uint8_t>& headerData, const Array<uint8_t>& data, const std::shared_ptr<Grid<Dimension, CoordinateType>>& loadedGrid, Serializer::DataToLoad dataToLoad) noexcept
		{
			AssertE(headerData.size() == sizeof(SerializedData::Header), MessageTag::Serialization, "Invalid deserialized problem header size!");
			const auto header = reinterpret_cast<typename SerializedData::Header*>(headerData.data());
			AssertE(header->dataType == static_cast<uint32_t>(Serializer::ToDataType<FieldType>()), MessageTag::Serialization, "Invalid deserialized problem data type!");
			AssertE(header->fieldCount == this->GetDescriptor().ContinuousEquationCount(), MessageTag::Serialization, "Invalid deserialized problem field count!");
			AssertE(header->variableCount == this->GetDescriptor().DiscreteEquationCount(), MessageTag::Serialization, "Invalid deserialized problem variable count!");
			AssertE(header->parameterCount == this->GetDescriptor().ParameterCount(), MessageTag::Serialization, "Invalid deserialized problem parameter count!");
			auto sourceData = reinterpret_cast<FieldType*>(data.data());
			if (dataToLoad & Serializer::DataToLoad::Parameters)
			{
				this->LoadParameters(sourceData);
			}
			sourceData += this->GetDescriptor().ParameterCount();
			if (dataToLoad & Serializer::DataToLoad::Time)
			{
				SetTime(*reinterpret_cast<CoordinateType*>(sourceData));
			}
			sourceData += sizeof(CoordinateType);
			if (dataToLoad & Serializer::DataToLoad::Variables)
			{
				this->LoadVariables(sourceData);
			}
		}
	};

	template<size_t Dimension, typename CoordinateType, typename FieldType>
	template<template<typename> typename MatrixType>
	[[nodiscard]] auto
		ExplicitTransientProblemDescriptor<Dimension, CoordinateType, FieldType>::MakeProblem
		( sptr<Grid<Dimension, CoordinateType>> grid
		, uptr<Discretization<Dimension, MatrixType, CoordinateType>> discretizer
		) const noexcept
	{
		AssertE(grid->GetDescriptor() == this->GetGridDescriptor() || this->Validate(),
			MessageTag::Problem, "Invalid problem descriptor.");
		return sptr<ExplicitTransientProblem<Dimension, MatrixType, CoordinateType, FieldType>>(
			new ExplicitTransientProblem<Dimension, MatrixType, CoordinateType, FieldType>(std::move(grid), std::move(discretizer), *this));
	}
}