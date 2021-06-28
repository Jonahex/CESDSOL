cmake_minimum_required(VERSION 3.10)

project(CESDSOL VERSION 1.0)

set(MainPath "" CACHE FILEPATH "Path to problem definition file containing main().")

set(MathBackend "MKLMath" CACHE STRING "Math backend to use.")
set_property(CACHE MathBackend PROPERTY STRINGS "MKLMath")

set(Parallelism "OpenMPParallelism" CACHE STRING "Parallelism backend to use.")
set_property(CACHE Parallelism PROPERTY STRINGS "SequentialParallelism")
set_property(CACHE Parallelism PROPERTY STRINGS "OpenMPParallelism")

set(DebugMode OFF CACHE BOOL "Add additional runtime checks.")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SourcesPath "CESDSOL")
file(GLOB_RECURSE Sources CONFIGURE_DEPENDS "${SourcesPath}/*.h" "${SourcesPath}/*.cpp")

foreach(Item IN ITEMS ${Sources})
    get_filename_component(ItemPath "${Item}" PATH)
    string(REPLACE "${CMAKE_SOURCE_DIR}" "" GroupPath "${ItemPath}")
    string(REPLACE "/" "\\" GroupPath "${GroupPath}")
    source_group("${GroupPath}" FILES "${Item}")
endforeach()

add_executable(CESDSOL ${Sources} ${MainPath})

find_package(OpenMP REQUIRED)
target_link_libraries(CESDSOL PRIVATE OpenMP::OpenMP_CXX)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	target_compile_options(CESDSOL PRIVATE "/openmp:llvm")
endif()

if (MathBackend STREQUAL "MKLMath")
	set(MKL_THREAD_LAYER "Intel OpenMP")
	find_package(MKL REQUIRED PATHS cmake/)
	target_link_libraries(CESDSOL PRIVATE MKL::Static)	
endif()

target_include_directories(CESDSOL PRIVATE ${SourcesPath})

target_compile_definitions(CESDSOL PRIVATE MathLibrary=${MathBackend} ParallelismBackend=${Parallelism} ExamplesPath="${CMAKE_CURRENT_LIST_DIR}/examples/")
if (DebugMode)
	target_compile_definitions(CESDSOL PRIVATE DebugMode)
endif()
