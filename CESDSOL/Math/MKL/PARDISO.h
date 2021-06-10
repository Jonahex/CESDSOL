#pragma once

#include "Math/MKL/Concepts.h"
#include "Math/MKL/CSRMatrix.h"
#include "Math/LinearAlgebra.h"
#include "Math/LinearSolver.h"
#include "Utils/Aliases.h"

#include "pardiso.h"

#include <chrono>
#include <concepts>

namespace CESDSOL::MKL
{
	template<MKLScalar ScalarType>
	class PARDISO final
		: public LinearSolver<CSRMatrix<ScalarType>, Array<ScalarType>>
	{
	private:
		using MatrixType = CSRMatrix<ScalarType>;
		static_assert(std::is_same_v<typename MatrixType::index_type, MKL_INT>);
		static_assert(MatrixType::StartingIndex == 0 || MatrixType::StartingIndex == 1);
		
		using VectorType = Vector<ScalarType>;

		enum class SolutionPhase
		{
			Analysis = 1, 
			AnalysisFactorization = 12, 
			AnalysisFactorizationSolveIterativeRefinement = 13,
			Factor = 22, 
			FactorizationSolveIterativeRefinement = 23,
			SolveIterativeRefinement = 33,
			SolveIterativeRefinementForward = 331,
			SolveIterativeRefinementDiagonal = 332,
			SolveIterativeRefinementBackward = 333,
			ReleaseLU = 0, 
			ReleaseAll = -1
		};

		enum class PardisoMatrixType 
		{ 
			RealStructurallySymmetric = 1, 
			RealSymmetricPositiveDefinite = 2, 
			RealSymmetricIndefinite = -2, 
			ComplexStructurallySymmetric = 3,
			ComplexHermiteanPositiveDefinite = 4,
			ComplexHermiteanIndefinite = -4, 
			ComplexSymmetric = 6, 
			RealNonsymmetric = 11, 
			ComplexNonsymmetric = 13 
		};

		void* internalData[64];
		MKL_INT intParameters[64];

		MKL_INT currentPhase;
		
		static constexpr MKL_INT EmptyInternalMatrix = -1;
		MKL_INT equationCount = EmptyInternalMatrix;
		
		bool useCgs = false;
		double cgsTolerance = 1e-6;
		Array<MKL_INT> permutation;

		static constexpr MKL_INT MklMatrixType = std::is_same_v<ScalarType, c32> || std::is_same_v<ScalarType, c64>
			? static_cast<MKL_INT>(PardisoMatrixType::ComplexNonsymmetric)
			: static_cast<MKL_INT>(PardisoMatrixType::RealNonsymmetric);
		static constexpr MKL_INT MaxFactorCount = 1;
		static constexpr MKL_INT MatrixNumber = 1;
		static constexpr MKL_INT RhsCount = 1;
		static constexpr MKL_INT MessageLevel = 0;

		void SetSolutionPhase(SolutionPhase phase) noexcept
		{
			currentPhase = static_cast<MKL_INT>(phase);
		}

		static const char* GetErrorMessage(MKL_INT error)
		{
			switch (error)
			{
			case -1:
				return "input inconsistent.";
			case -2:
				return "not enough memory.";
			case -3:
				return "reordering problem.";
			case -4:
				return "zero pivot, numerical factorization or iterative refinement problem.";
			case -6:
				return "reordering failed.";
			case -7:
				return "diagonal matrix is singular.";
			case -8:
				return "32-bit integer overflow problem.";
			case -9:
				return "not enough memory for OOC.";
			case -10:
				return "problems with opening OOC temporary files.";
			case -11:
				return "read/write problems with the OOC data file.";
			case -12:
				return "pardiso_64 called from 32-bit library.";
			case -13:
				return "interrupted by the mkl_progress function.";
			default:
				return "unclassified (internal) error.";
			}
		}

		static void NotifyError(MKL_INT error) noexcept
		{
			if (error < 0)
			{
				Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver, 
					Format("PARDISO encountered error: {}", GetErrorMessage(error)));
			}
		}

		void UpdateCgsOptions() noexcept
		{
			MKL_INT value = -10 * std::round(std::log10(cgsTolerance));
			if (value < 0)
			{
				Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::LinearSolver, 
					Format("Trying to set invalid CGS tolerance value {} for PARDISO solver.", cgsTolerance));
			}
			else
			{
				if (useCgs)
				{
					value += 1;
				}
				intParameters[3] = value;
			}
		}

	public:
		enum class FillInReductionMethod
		{
			MinimumDegree = 0,
			NestedDissection = 2,
			ParallelNestedDissection = 3
		};

		enum class PermutationMode
		{
			Ignore = 0,
			UserProvided = 1,
			Calculated = 2
		};

		enum class FactorizationMethod
		{
			Classic = 0,
			TwoLevel = 1,
			ImprovedTwoLevel = 10
		};

		enum class SolveParallelizationStrategy
		{
			ParallelOnRHSVectors = 0,
			Sequential = 1,
			Parallel = 2
		};

		void SetFillInReductionMethod(FillInReductionMethod method) noexcept
		{
			intParameters[1] = static_cast<MKL_INT>(method);
		}

		void SetCGSTolerance(double value) noexcept
		{
			cgsTolerance = value;
			UpdateCgsOptions();
		}

		void SetUseCGS(bool value) noexcept
		{
			useCgs = value;
			UpdateCgsOptions();
		}

		template<typename ArrayReference> requires (std::same_as<std::remove_cvref_t<ArrayReference>, Array<MKL_INT>>)
		void SetUserPermutation(ArrayReference&& array) noexcept
		{
			permutation = std::forward<ArrayReference>(array);
		}

		void SetPermutationMode(PermutationMode mode) noexcept
		{
			intParameters[4] = static_cast<MKL_INT>(mode);
		}

		[[nodiscard]] const Array<MKL_INT>& GetPermutation() const noexcept
		{
			return permutation;
		}

		void SetPivotPerturbationFactor(double value) noexcept
		{
			if (value > 1 || value <= 0)
			{
				Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::LinearSolver, Format("Trying to set invalid pivot perturbation factor {} for PARDISO solver.", value));
			}
			else
			{
				intParameters[9] = std::round(-std::log10(value));
			}
		}

		void SetUseDiagonalElementScaling(bool ifUse) noexcept
		{
			if (ifUse && GetFactorizationMethod() != FactorizationMethod::Classic)
			{
				Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::LinearSolver,
					"Diagonal scaling is incompatible with two-level factorization!");
			}
			else
			{
				intParameters[10] = static_cast<MKL_INT>(ifUse);
			}
		}

		[[nodiscard]] bool UseDiagonalElementScaling() const noexcept
		{
			return static_cast<bool>(intParameters[10]);
		}

		void SetUseWeightedMatching(bool ifUse) noexcept
		{
			if (ifUse && GetFactorizationMethod() != FactorizationMethod::Classic)
			{
				Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::LinearSolver,
					"Weighted matching is incompatible with two-level factorization!");
			}
			else
			{
				intParameters[12] = static_cast<MKL_INT>(ifUse);
			}
		}

		[[nodiscard]] bool UseWeightedMatching() const noexcept
		{
			return static_cast<bool>(intParameters[10]);
		}

		void SetFactorizationMethod(FactorizationMethod method) noexcept
		{
			if (method != FactorizationMethod::Classic && (UseDiagonalElementScaling() || UseWeightedMatching()))
			{
				Logger::Log(MessageType::Warning, MessagePriority::High, MessageTag::LinearSolver, 
					"Two-level factorization requires both diagonal element scaling and weighted matching disabled!");
			}
			else
			{
				intParameters[23] = static_cast<MKL_INT>(method);
			}
		}

		[[nodiscard]] FactorizationMethod GetFactorizationMethod() const noexcept
		{
			return static_cast<FactorizationMethod>(intParameters[23]);
		}

		void SetSolveParallelizationStrategy(SolveParallelizationStrategy strategy) noexcept
		{
			intParameters[24] = static_cast<MKL_INT>(strategy);
		}

		PARDISO() noexcept
		{
			pardisoinit(internalData, &MklMatrixType, intParameters);

			SetFillInReductionMethod(FillInReductionMethod::ParallelNestedDissection);
			SetSolveParallelizationStrategy(SolveParallelizationStrategy::Parallel);
			SetUseCGS(true);
			if constexpr (std::is_same_v<ScalarType, f32> || std::is_same_v<ScalarType, c32>)
			{
				intParameters[27] = 1;
			}
			intParameters[34] = static_cast<MKL_INT>(MatrixType::StartingIndex == 0);

			SetSolutionPhase(SolutionPhase::AnalysisFactorizationSolveIterativeRefinement);
		};

		~PARDISO() noexcept
		{
			ResetSolutionData();
		}

		bool Solve(const MatrixType& matrix, const VectorType& y, VectorType& x) noexcept override
		{
			Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::LinearSolver, Format("Starting solving system of {} linear equations with PARDISO.", matrix.RowCount()));
			std::chrono::high_resolution_clock clock;
			const auto solutionStartTime = clock.now();
			MKL_INT error;
			if (equationCount != EmptyInternalMatrix)
			{
				AssertE(equationCount == matrix.RowCount(), MessageTag::LinearSolver, "Inconsistent nonzero count in matrix for PARDISO to solve.");
			}
			else
			{
				equationCount = static_cast<MKL_INT>(matrix.RowCount());
			}
			pardiso(internalData, &MaxFactorCount, &MatrixNumber, &MklMatrixType, &currentPhase, &equationCount,
				matrix.GetValues().data(), matrix.GetRowCounts().data(), matrix.GetColumnIndices().data(),
				permutation.data(), &RhsCount, intParameters, &MessageLevel,
				y.data(), x.data(), &error);
			NotifyError(error);
			if (error == 0)
			{
				SetSolutionPhase(SolutionPhase::FactorizationSolveIterativeRefinement);
				Logger::Log(MessageType::Info, MessagePriority::Medium, MessageTag::LinearSolver, 
					Format("Linear system is solved by PARDISO in {}", clock.now() - solutionStartTime).c_str());
			}
			return error == 0;
		};

		void ResetSolutionData() noexcept
		{
			if (equationCount != EmptyInternalMatrix)
			{
				MKL_INT error;
				SetSolutionPhase(SolutionPhase::ReleaseAll);
				pardiso(internalData, &MaxFactorCount, &MatrixNumber, &MklMatrixType, &currentPhase, &equationCount, nullptr,
					nullptr, nullptr, permutation.data(), &RhsCount, intParameters, &MessageLevel,
					nullptr, nullptr, &error);
				NotifyError(error);
				SetSolutionPhase(SolutionPhase::AnalysisFactorizationSolveIterativeRefinement);
				equationCount = EmptyInternalMatrix;
			}
		}
	};
}