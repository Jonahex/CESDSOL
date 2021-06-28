#pragma once

#include "Configuration.h"

#define SingleArg(...) __VA_ARGS__

#if ParallelismBackend == OpenMPParallelism
#include "Utils/Parallelism/OpenMP.h"
#elif ParallelsimBackend == SequentialParallelism
#include "Utils/Parallelism/Sequential.h"
#endif 