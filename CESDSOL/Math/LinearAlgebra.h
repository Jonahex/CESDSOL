#pragma once

#include "../Configuration.h"

#if MathLibrary == NativeMath
#include "Array.h"
#include "DenseMatrix.h"
#include "CSRMatrix.h"
#include "Native/VectorOperations.h"

namespace CESDSOL
{
	template<typename ScalarType>
	using Vector = Array<ScalarType>;
	template<typename ScalarType>
	using CSRMatrix = Native::CSRMatrix<ScalarType, size_t, 0>;

	namespace LinearAlgebra = Native;
}
#elif MathLibrary == MKLMath
#include "Array.h"
#include "DenseMatrix.h"
#include "MKL/CSRMatrix.h"
#include "MKL/CSRMatrixOperations.h"
#include "MKL/VectorOperations.h"

namespace CESDSOL
{
	template<typename ScalarType>
	using Vector = Array<ScalarType>;
	template<typename ScalarType>
	using CSRMatrix = MKL::CSRMatrix<ScalarType>;

	namespace LinearAlgebra = MKL;
}
#endif

#include "CSRMatrixOperations.h"
#include "SparseVector.h" 
#include "SparseVectorOperations.h"