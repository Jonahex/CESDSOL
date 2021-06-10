#pragma once

#include "Math/MKL/CSRMatrix.h"
#include "Math/MKL/Utils.h"
#include "Math/Concepts.h"

namespace CESDSOL
{
	template<MKL::MKLCSRMatrix MatrixType, Concepts::Vector VectorType, typename ScalarType = f64>
	void Multiply(const MatrixType& A, const VectorType& x, VectorType& y, ScalarType alpha = 1., ScalarType beta = 0.)
	{
		AssertE(A.ColumnCount() == x.size(), MessageTag::Math, "Trying to multiply matrix and vector with incompatible sizes.");
		AssertE(x.size() == y.size(), MessageTag::Math, "Trying to assign vectors with incompatible sizes.");

#define SPARSE_MV_OPERATION(x) mkl_sparse_##x##_mv
		MKL_SPARSE_OPERATION(SPARSE_MV_OPERATION, "Sparse matrix-vector multiplication", SPARSE_OPERATION_NON_TRANSPOSE,
			alpha, A.GetHandle(), A.GetDescriptor(), x.data(), beta, y.data());
#undef SPARSE_MV_OPERATION
	}

	template<MKL::MKLCSRMatrix MatrixType, Concepts::Vector VectorType, typename ScalarType = f64>
	VectorType Multiply(const MatrixType& A, const VectorType& x, ScalarType alpha = 1.)
	{
		AssertE(A.ColumnCount() == x.size(), MessageTag::Math, "Trying to multiply matrix and vector with incompatible sizes.");

		auto y = VectorType(A.RowCount());
		Multiply(A, x, y, alpha);
		return y;
	}

	template<MKL::MKLCSRMatrix MatrixType>
	MatrixType Multiply(const MatrixType& A, const MatrixType& B)
	{
		AssertE(A.ColumnCount() == B.RowCount(), MessageTag::Math, "Trying to multiply matrices with incompatible sizes.");

		sparse_matrix_t newHandle;
		auto status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A.GetHandle(), B.GetHandle(), &newHandle);
		AssertE(status == SPARSE_STATUS_SUCCESS, MessageTag::Math, Format("Sparse matrix-matrix multiplication: {}.", MKL::SparseStatusMessages[status]));

		return MatrixType(newHandle);
	}

	template<MKL::MKLCSRMatrix MatrixType, typename ScalarType = f64>
	MatrixType Add(const MatrixType& A, const MatrixType& B, ScalarType alpha = 1.)
	{
		AssertE(A.RowCount() == B.RowCount() && A.ColumnCount() == B.ColumnCount(), MessageTag::Math,
			"Trying to add matrices with incompatible sizes.");

		sparse_matrix_t newHandle;
#define SPARSE_ADD_OPERATION(x) mkl_sparse_##x##_add
		MKL_SPARSE_OPERATION(SPARSE_ADD_OPERATION, "Sparse matrix addition", SPARSE_OPERATION_NON_TRANSPOSE,
			A.GetHandle(), alpha, B.GetHandle(), &newHandle);
#undef SPARSE_ADD_OPERATION
		return MatrixType(newHandle);
	}

	template<MKL::MKLCSRMatrix MatrixType, Concepts::Vector XVectorType, Concepts::Vector YVectorType, typename ScalarType = f64>
	void TriangularSolve(const MatrixType& A, const XVectorType& x, YVectorType& y, ScalarType alpha = 1., bool isUpper = true)
	{
		AssertE(A.ColumnCount() == x.size() && A.RowCount() == y.size(), MessageTag::Math, 
			"Trying to perform triangular solve for matrix and vector with incompatible sizes.");
		
		auto descriptor = A.GetDescriptor();
		descriptor.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
		if (isUpper)
		{
			descriptor.mode = SPARSE_FILL_MODE_UPPER;
			descriptor.diag = SPARSE_DIAG_NON_UNIT;
		}
		else
		{
			descriptor.mode = SPARSE_FILL_MODE_LOWER;
			descriptor.diag = SPARSE_DIAG_UNIT;
		}
		
#define SPARSE_TRSV_OPERATION(x) mkl_sparse_##x##_trsv
		MKL_SPARSE_OPERATION(SPARSE_TRSV_OPERATION, "Sparse matrix triangular solve", SPARSE_OPERATION_NON_TRANSPOSE,
			alpha, A.GetHandle(), descriptor, x.data(), y.data());
#undef SPARSE_TRSV_OPERATION
	}
}