#pragma once

#include "Math/Preconditioner.h"

#include "krylov/krylov.h"

namespace CESDSOL::HYPRE
{
	template<typename MatrixType, typename ArrayType>
	class BiCGSTAB final
		: public LinearSolver<MatrixType, ArrayType>
	{
	public:
		BiCGSTAB(uptr<Preconditioner<MatrixType, ArrayType>> aPreconditioner = nullptr)
			: preconditioner(std::move(aPreconditioner))
		{
			hypre_BiCGSTABFunctions* solverFunctions = hypre_BiCGSTABFunctionsCreate(&CreateVector, &DestroyVector, &MatvecCreate, &Matvec, &MatvecDestroy,
				&InnerProd, &CopyVector, &ClearVector, &ScaleVector, &Axpy, &CommInfo,
				&PrecondSetup, &Precond);
			solver = hypre_BiCGSTABCreate(solverFunctions);
			hypre_BiCGSTABSetLogging(solver, 1);
			hypre_BiCGSTABSetPrintLevel(solver, 1);
		}

		~BiCGSTAB()
		{
			hypre_BiCGSTABDestroy(solver);
		}

		bool Solve(const MatrixType& matrix, const ArrayType& y, ArrayType& x) override
		{
			hypre_BiCGSTABSetPrecond(solver, &Precond, &PrecondSetup, preconditioner.get());
			hypre_BiCGSTABSetup(solver, const_cast<MatrixType*>(&matrix), const_cast<ArrayType*>(&y), &x);
			hypre_BiCGSTABSolve(solver, const_cast<MatrixType*>(&matrix), const_cast<ArrayType*>(&y), &x);
			return true;
		}

		void SetPreconditioner(uptr<Preconditioner<MatrixType, ArrayType>> aPreconditioner) noexcept
		{
			preconditioner = std::move(aPreconditioner);
		}

	private:
		static void* CreateVector(void* vvector)
		{
			return new ArrayType(static_cast<ArrayType*>(vvector)->size());
		}

		static HYPRE_Int DestroyVector(void* vvector)
		{
			delete vvector;
			return 0;
		}

		static void* MatvecCreate(void* A, void* x)
		{
			return nullptr;
		}

		static HYPRE_Int Matvec(void* matvec_data, HYPRE_Complex alpha, void* A,
			void* x, HYPRE_Complex beta, void* y)
		{
			MVMultiply(*static_cast<MatrixType*>(A), *static_cast<ArrayType*>(x), *static_cast<ArrayType*>(y), alpha, beta);
			return 0;
		}

		static HYPRE_Int MatvecDestroy(void* matvec_data)
		{
			return 0;
		}

		static HYPRE_Real InnerProd(void* x, void* y)
		{
			return DotProduct(*static_cast<ArrayType*>(x), *static_cast<ArrayType*>(y));
		}

		static HYPRE_Int CopyVector(void* x, void* y)
		{
			Copy(*static_cast<ArrayType*>(x), *static_cast<ArrayType*>(y));
			return 0;
		}

		static HYPRE_Int ClearVector(void* x)
		{
			Fill(*static_cast<ArrayType*>(x), 0.);
			return 0;
		}

		static HYPRE_Int ScaleVector(HYPRE_Complex alpha, void* x)
		{
			Scale(alpha, *static_cast<ArrayType*>(x));
			return 0;
		}

		static HYPRE_Int Axpy(HYPRE_Complex alpha, void* x, void* y)
		{
			AXPY(alpha, *static_cast<ArrayType*>(x), *static_cast<ArrayType*>(y));
			return 0;
		}

		static HYPRE_Int CommInfo(void* A, HYPRE_Int* my_id,
			HYPRE_Int* num_procs)
		{
			*my_id = 0;
			*num_procs = 1;
			return 0;
		}

		static HYPRE_Int PrecondSetup(void* vdata, void* A, void* b, void* x)
		{
			if (vdata != nullptr)
			{
				static_cast<Preconditioner<MatrixType, ArrayType>*>(vdata)->Setup(*static_cast<MatrixType*>(A), *static_cast<ArrayType*>(b));
			}
			return 0;
		}

		static HYPRE_Int Precond(void* vdata, void* A, void* b, void* x)
		{
			if (vdata != nullptr)
			{
				static_cast<Preconditioner<MatrixType, ArrayType>*>(vdata)->Solve(*static_cast<MatrixType*>(A), *static_cast<ArrayType*>(b), *static_cast<ArrayType*>(x));
			}
			else
			{
				Copy(*static_cast<ArrayType*>(b), *static_cast<ArrayType*>(x));
			}
			return 0;
		}

		uptr<Preconditioner<MatrixType, ArrayType>> preconditioner;

		void* solver;
	};
}