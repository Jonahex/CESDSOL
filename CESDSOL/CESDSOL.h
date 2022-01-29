#pragma once

#include "Configuration.h"
#include "Discretization/StructuredFiniteDifferenceDiscretization.h"
#include "Grid/DirectProductGrid.h"
#include "Grid/Grid.h"
#include "Math/GoldenSectionSearch.h"
#include "Math/ModifiedNewton.h"
#include "Math/ODE/Tables/BogackiShampine32.h"
#include "Math/ODE/Tables/DormandPrince54.h"
#include "Math/ODE/Tables/DormandPrince853.h"
#include "Math/ODE/Tables/Euler1.h"
#include "Math/ODE/Tables/Feagin109.h"
#include "Math/ODE/Tables/Fehlberg21.h"
#include "Math/ODE/Tables/Heun21.h"
#include "Math/ODE/Tables/Midpoint21.h"
#include "Math/ODE/Tables/OwrenZennaro32.h"
#include "Math/ODE/Tables/OwrenZennaro43.h"
#include "Math/ODE/Tables/OwrenZennaro54.h"
#include "Math/ODE/Tables/Ralston21.h"
#include "Math/ODE/Tables/RungeKutta4.h"
#include "Math/ODE/Tables/TanakaYamashita76.h"
#include "Math/ODE/Tables/Tsitouras54.h"
#include "Math/ODE/Tables/TsitourasPapakostas87.h"
#include "Math/ODE/Tables/Verner87.h"
#include "Math/ODE/RungeKuttaSolver.h"
#include "Math/TrivialLineSearcher.h"
#include "Math/VectorOperations.h"
#include "ParametricSweep/AdaptiveParametricSweeper.h"
#include "ParametricSweep/FixedStepParametricSweeper.h"
#include "Problem/ExplicitTransientProblem.h"
#include "Problem/StationaryProblem.h"
#include "Serialization/Serializer.h"

#if MathLibrary == MKLMath
#include "Math/MKL/FGMRES.h"
#include "Math/MKL/ILU0.h"
#include "Math/MKL/ILUT.h"
#include "Math/MKL/PARDISO.h"
#endif

#ifdef UseHYPRE
#include "Math/HYPRE/BiCGSTAB.h"
#endif