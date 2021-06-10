#pragma once

#include "Math/MKL/Concepts.h"
#include "Math/MKL/Utils.h"
#include "Math/CSRMatrix.h"

#include "mkl.h"

#include <concepts>

namespace CESDSOL::MKL
{
	template<MKLScalar ScalarType, MKL_INT StartingIndex = 1>
	class CSRMatrix final
		: public Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>
	{
	private:
		sparse_matrix_t handle = nullptr;
		matrix_descr descriptor;

		void ExportHandle() noexcept
		{
			ScalarType* values;
			MKL_INT* columnIndices, * rowCounts, * rowEnds;
			MKL_INT rowCount, columnCount;
			sparse_index_base_t startingIndex;
#define EXPORT_CSR_OPERATION(x) mkl_sparse_##x##_export_csr
			MKL_SPARSE_OPERATION(EXPORT_CSR_OPERATION, "CSR matrix export", handle, &startingIndex, &rowCount, &columnCount,
				&rowCounts, &rowEnds, &columnIndices, &values);
#undef EXPORT_CSR_OPERATION
			AssertE(startingIndex == StartingIndex, MessageTag::Math, "Incompatible starting indices of CSR matrix");
			this->rowCount = rowCount;
			this->columnCount = columnCount;
			this->rowCounts = Array<MKL_INT>(rowCounts, this->RowCount() + 1);
			this->columnIndices = Array<MKL_INT>(columnIndices, this->NonZeroCount());
			this->values = Array<ScalarType>(values, this->NonZeroCount());
		}

		void CreateHandle() noexcept
		{
#define CREATE_CSR_OPERATION(x) mkl_sparse_##x##_create_csr
			MKL_SPARSE_OPERATION(CREATE_CSR_OPERATION, "CSR matrix construction",
				&handle, static_cast<sparse_index_base_t>(StartingIndex),
				static_cast<int>(this->rowCount), static_cast<int>(this->columnCount),
				this->rowCounts.data(), this->rowCounts.data() + 1, this->columnIndices.data(), this->values.data());
#undef CREATE_CSR_OPERATION
		}

	public:
		constexpr CSRMatrix() noexcept = default;

		constexpr CSRMatrix(MKL_INT aRowCount, MKL_INT aColumnCount, MKL_INT nonZeroCount) noexcept
			: Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>(aRowCount, aColumnCount, nonZeroCount)
			, descriptor({ SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_FULL, SPARSE_DIAG_NON_UNIT })
		{
			CreateHandle();
		}

		constexpr CSRMatrix(sparse_matrix_t aHandle) noexcept
			: handle(aHandle)
			, descriptor({ SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_FULL, SPARSE_DIAG_NON_UNIT })
		{
			ExportHandle();
		}

		constexpr CSRMatrix(CSRMatrix<ScalarType, StartingIndex>&& other) noexcept
			: Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>(std::move(other))
			, handle(other.handle)
			, descriptor(other.descriptor)
		{
			other.handle = nullptr;
		}

		constexpr CSRMatrix(const CSRMatrix<ScalarType, StartingIndex>& other) noexcept
			: Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>(other)
			, descriptor(other.descriptor)
		{
			CreateHandle();
		}

		constexpr CSRMatrix<ScalarType, StartingIndex>& operator=(CSRMatrix<ScalarType, StartingIndex>&& other) noexcept
		{
			if (this != &other)
			{
				Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>::operator=(std::move(other));
				handle = other.handle;
				descriptor = other.descriptor;
				other.handle = nullptr;
			}
			return *this;
		}

		constexpr CSRMatrix<ScalarType, StartingIndex>& operator=(const CSRMatrix<ScalarType, StartingIndex>& other) noexcept
		{
			if (this != &other)
			{
				Native::CSRMatrix<ScalarType, MKL_INT, StartingIndex>::operator=(other);
				descriptor = other.descriptor;
				CreateHandle();
			}
			return *this;
		}

		[[nodiscard]] sparse_matrix_t GetHandle() const noexcept
		{
			return handle;
		}

		[[nodiscard]] const matrix_descr& GetDescriptor() const noexcept
		{
			return descriptor;
		}

		[[nodiscard]] matrix_descr& GetDescriptor() noexcept
		{
			return descriptor;
		}

		~CSRMatrix()
		{
			if (handle != nullptr)
			{
				this->values.Release();
				this->rowCounts.Release();
				this->columnIndices.Release();
				
				const auto status = mkl_sparse_destroy(handle);
				AssertE(status == SPARSE_STATUS_SUCCESS, MessageTag::Math, Format("Sparse matrix destruction: {}.", MKL::SparseStatusMessages[status]));
			}
		}
	};

	template<typename MatrixType>
	concept MKLCSRMatrix = std::same_as<MatrixType, CSRMatrix<typename MatrixType::value_type, MatrixType::StartingIndex>>;
}