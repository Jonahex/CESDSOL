#pragma once

#include "Problem/StationaryProblem.h"

namespace CESDSOL::SimplifiedInterface
{
	template<typename LocalsType>
	struct FieldHelper
	{
		const LocalsType& locals;
		size_t fieldIndex;

		FieldHelper(const LocalsType& aLocals, size_t aFieldIndex)
			: locals(aLocals)
			, fieldIndex(aFieldIndex)
		{}

		constexpr operator typename LocalsType::FieldType() const
		{
			return locals.FieldValues[fieldIndex];
		}
	};

	template<typename LocalsType, typename GlobalsType>
	struct VariableHelper
	{
		const GlobalsType& globals;
		size_t variableIndex;
		size_t fieldIndex;

		VariableHelper(const LocalsType& locals, const GlobalsType& aGlobals, size_t aVariableIndex)
			: globals(aGlobals)
			, variableIndex(aVariableIndex)
			, fieldIndex(locals.FieldValues.size() + variableIndex)
		{}

		constexpr operator typename GlobalsType::FieldType() const
		{
			return globals.DiscreteVariables[variableIndex];
		}
	};

	template<typename LocalsType>
	struct DirectionHelper
	{
		const LocalsType& locals;
		size_t directionIndex;

		DirectionHelper(const LocalsType& aLocals, size_t aDirectionIndex)
			: locals(aLocals)
			, directionIndex(aDirectionIndex)
		{}

		constexpr operator typename LocalsType::FieldType() const
		{
			return locals.Point[directionIndex];
		}
	};

	template<typename LocalsType>
	struct DerivativeHelper
	{
		const LocalsType& locals;
		size_t fieldIndex;
		size_t operatorIndex;

		DerivativeHelper(const LocalsType& aLocals, const FieldHelper<LocalsType>& fieldHelper, size_t aOperatorIndex)
			: locals(aLocals)
			, fieldIndex(fieldHelper.fieldIndex)
			, operatorIndex(aOperatorIndex)
		{}

		constexpr operator typename LocalsType::FieldType() const
		{
			return locals.DerivativeValues[fieldIndex][operatorIndex];
		}
	};

	struct Void {};

	template<>
	struct FieldHelper<Void>
	{
		size_t fieldIndex;
		static constexpr int operatorIndex = -1;

		FieldHelper(Void aLocals, size_t aFieldIndex)
			: fieldIndex(aFieldIndex)
		{}
	};

	template<>
	struct VariableHelper<Void, Void>
	{
		size_t fieldIndex;
		static constexpr int operatorIndex = -1;

		VariableHelper(Void aLocals, Void aGlobals, size_t aVariableIndex)
			: fieldIndex(aVariableIndex)
		{}
	};

	template<>
	struct DirectionHelper<Void>
	{
		size_t directionIndex;

		DirectionHelper(Void aLocals, size_t aDirectionIndex)
			: directionIndex(aDirectionIndex)
		{}
	};

	template<>
	struct DerivativeHelper<Void>
	{
		size_t fieldIndex;
		size_t operatorIndex;

		DerivativeHelper(Void aLocals, const FieldHelper<Void>& fieldHelper, size_t aOperatorIndex)
			: fieldIndex(fieldHelper.fieldIndex)
			, operatorIndex(aOperatorIndex)
		{}
	};
}

inline constexpr CESDSOL::SimplifiedInterface::Void SimplifiedInterfaceLocalsHelper;
inline constexpr CESDSOL::SimplifiedInterface::Void SimplifiedInterfaceGlobalsHelper;

#define Field(index) CESDSOL::SimplifiedInterface::FieldHelper(SimplifiedInterfaceLocalsHelper, index)
#define DerivativeOperator(field, operatorIndex) CESDSOL::SimplifiedInterface::DerivativeHelper(SimplifiedInterfaceLocalsHelper, field, operatorIndex)
#define Equation(index) index
#define Region(index) index
#define Interior Region(0)
#define Direction(index) CESDSOL::SimplifiedInterface::DirectionHelper(SimplifiedInterfaceLocalsHelper, index)
#define LeftBoundary(index) Region(2 * index.directionIndex + 1)
#define RightBoundary(index) Region(2 * index.directionIndex + 2)
#define AddFieldEquation(problem, equation, region, rhs) problem.SetContinuousEquation(equation, region, [](const auto& SimplifiedInterfaceLocalsHelper, const auto& SimplifiedInterfaceGlobalsHelper) {return rhs;})
#define AddVariableEquation(problem, equation, rhs) problem.SetDiscreteEquation(equation, [](const auto& SimplifiedInterfaceGlobalsHelper) {return rhs;})
#define SetJacobian(problem, equation, operator, region, rhs) problem.SetJacobianComponent(equation, operator.fieldIndex, operator.operatorIndex + 1, region, [](const auto& SimplifiedInterfaceLocalsHelper, const auto& SimplifiedInterfaceGlobalsHelper) {return rhs;})
#define Parameter(index) SimplifiedInterfaceGlobalsHelper.Parameters[index]
#define Variable(index) CESDSOL::SimplifiedInterface::VariableHelper(SimplifiedInterfaceLocalsHelper, SimplifiedInterfaceGlobalsHelper, index)