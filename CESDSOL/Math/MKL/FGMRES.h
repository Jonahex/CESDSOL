#pragma once

#include "Math/LinearAlgebra.h"
#include "Math/LinearSolver.h"
#include "Math/Preconditioner.h"

namespace CESDSOL::MKL
{
	enum class FGMRESExitConditions : uint32_t
	{
		ResidualGoalReached = 1 << 0,
		IterationCount = 1 << 1,
		ZeroSolutionNorm = 1 << 2,
		Everything = (1 << 3) - 1
	};
	MakeFlag(FGMRESExitConditions)
	
	template<typename MatrixType>
	class FGMRES final
		: public LinearSolver<MatrixType, Vector<double>>
	{
	public:
		FGMRES(uptr<Preconditioner<MatrixType>> aPreconditioner = nullptr)
			: preconditioner(std::move(aPreconditioner))
		{}
		
		bool Solve(const MatrixType& matrix, const Vector<double>& y, Vector<double>& x) override
		{
			IntParametersType intParameters;
			DoubleParametersType floatParameters;

			MKL_INT rciRequest = 0;
			const MKL_INT size = y.size();
			auto tmp = MakeTempVector(size);
			dfgmres_init(&size, x.data(), y.data(), &rciRequest, intParameters.data(), 
				floatParameters.data(), tmp.data());
			InitParameters(intParameters, floatParameters);

			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LinearSolver, 
				Format("Starting solving system of {} linear equations with FGMRES.", size));
			dfgmres_check(&size, x.data(), y.data(), &rciRequest, intParameters.data(),
				floatParameters.data(), tmp.data());
			if (rciRequest != SuccessRCIRequest)
			{
				Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
					"Invalid parameters were passed to FGMRES.");
				return false;
			}

			Vector<double> preconditionerTmp;
			MatrixType preconditionedMatrix;
			if (preconditioner != nullptr)
			{
				preconditionerTmp = Vector<double>(size);
				preconditionedMatrix = matrix;
				if (!preconditioner->Apply(preconditionedMatrix)->success)
				{
					Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
						"Failed to calculate preconditioner matrix for FGMRES.");
					return false;
				}
			}
			while (true)
			{
				dfgmres(&size, x.data(), y.data(), &rciRequest, intParameters.data(),
					floatParameters.data(), tmp.data());
				if (rciRequest == MultiplyMatrixRCIRequest)
				{
					const std::span inSpan(tmp.data() + intParameters[21] - 1, size);
					std::span outSpan(tmp.data() + intParameters[22] - 1, size);
					Multiply(matrix, inSpan, outSpan, 1., 0.);
					continue;
				}
				if (rciRequest == ApplyPreconditionerRCIRequest)
				{
					const std::span inSpan(tmp.data() + intParameters[21] - 1, size);
					std::span outSpan(tmp.data() + intParameters[22] - 1, size);
					TriangularSolve(preconditionedMatrix, inSpan, preconditionerTmp, 1., false);
					TriangularSolve(preconditionedMatrix, preconditionerTmp, outSpan, 1., true);
					continue;
				}
				if (rciRequest == SuccessRCIRequest)
				{
					break;
				}
				if (rciRequest == HitIterationLimitWithoutConvergenceRCIRequest)
				{
					Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LinearSolver,
						"FGMRES reached iteration limit, but relative tolerance was not satisfied!");
					break;
				}
				if (rciRequest == DividedByZeroRCIRequest)
				{
					Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
						"FGMRES stopped due to division by zero!");
					return false;
				}
				if (rciRequest == InfiniteCycleRCIRequest)
				{
					Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
						"FGMRES stopped due to entering infinite cycle!");
					return false;
				}
				if (rciRequest == ErrorInParametersRCIRequest)
				{
					Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
						"FGMRES stopped due to parameters error!");
					return false;
				}
				if (rciRequest != CheckNormRCIRequest)
				{
					Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
						"FGMRES internal error!");
					return false;
				}
			}
			int iterationCount;
			dfgmres_get(&size, x.data(), y.data(), &rciRequest, intParameters.data(),
				floatParameters.data(), tmp.data(), &iterationCount);
			if (rciRequest != SuccessRCIRequest)
			{
				Logger::Log(MessageType::Error, MessagePriority::High, MessageTag::LinearSolver,
					"FGMRES internal error on data retrieval!");
				return false;
			}
			Logger::Log(MessageType::Info, MessagePriority::High, MessageTag::LinearSolver,
				Format("FGMRES solved linear system in {} iterations.", iterationCount));
			return true;
		}

		void SetPreconditioner(uptr<Preconditioner<MatrixType>> aPreconditioner) noexcept
		{
			preconditioner = std::move(aPreconditioner);
		}

	private:
		static constexpr MKL_INT ErrorInParametersRCIRequest = -12;
		static constexpr MKL_INT InfiniteCycleRCIRequest = -11;
		static constexpr MKL_INT DividedByZeroRCIRequest = -10;
		static constexpr MKL_INT HitIterationLimitWithoutConvergenceRCIRequest = -1;
		static constexpr MKL_INT SuccessRCIRequest = 0;
		static constexpr MKL_INT MultiplyMatrixRCIRequest = 1;
		static constexpr MKL_INT StoppingTestRCIRequest = 2;
		static constexpr MKL_INT ApplyPreconditionerRCIRequest = 3;
		static constexpr MKL_INT CheckNormRCIRequest = 4;

		using IntParametersType = std::array<MKL_INT, 128>;
		using DoubleParametersType = std::array<double, 128>;

		[[nodiscard]] Vector<double> MakeTempVector(MKL_INT problemSize) const noexcept
		{
			return Vector<double>((2 * restartIterationLimit + 1) * problemSize + restartIterationLimit * (restartIterationLimit + 9) / 2 + 1);
		}

		void InitParameters(IntParametersType& intParameters, DoubleParametersType& doubleParameters) const noexcept
		{
			intParameters[4] = iterationLimit;
			intParameters[5] = 0;
			intParameters[6] = 0;
			intParameters[7] = exitConditions & FGMRESExitConditions::IterationCount;
			intParameters[8] = exitConditions & FGMRESExitConditions::ResidualGoalReached;
			intParameters[9] = 0;
			intParameters[10] = preconditioner != nullptr;
			intParameters[11] = exitConditions & FGMRESExitConditions::ZeroSolutionNorm;
			intParameters[14] = restartIterationLimit;

			doubleParameters[0] = relativeTolerance;
			doubleParameters[1] = absoluteTolerance;
			doubleParameters[7] = zeroNormTolerance;
		}

		uptr<Preconditioner<MatrixType>> preconditioner;

		MakeProperty(exitConditions, ExitConditions, FGMRESExitConditions, FGMRESExitConditions::Everything)
		MakeProperty(iterationLimit, IterationLimit, size_t, 150)
		MakeProperty(restartIterationLimit, RestartIterationLimit, size_t, 150)
		MakeProperty(relativeTolerance, RelativeTolerance, double, 1e-6)
		MakeProperty(absoluteTolerance, AbsoluteTolerance, double, 0)
		MakeProperty(zeroNormTolerance, ZeroNormTolerance, double, 1e-12)
	};
}