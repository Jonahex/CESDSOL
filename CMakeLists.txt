cmake_minimum_required(VERSION 3.10)

project(CESDSOL VERSION 1.0)

set(MainPath "" CACHE FILEPATH "Path to problem definition file containing main().")

set(MathBackend "MKLMath" CACHE STRING "Math backend to use.")
set_property(CACHE MathBackend PROPERTY STRINGS "MKLMath")

set(Parallelism "OpenMPParallelism" CACHE STRING "Parallelism backend to use.")
set_property(CACHE Parallelism PROPERTY STRINGS "SequentialParallelism")
set_property(CACHE Parallelism PROPERTY STRINGS "OpenMPParallelism")

option(UseHYPRE "Use HYPRE library." OFF)
option(DebugMode "Add additional runtime checks." OFF)

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
	set(MKL_LINK "static")
	set(MKL_INTERFACE "lp64")
	find_package(MKL)	
	target_link_libraries(CESDSOL PRIVATE MKL::MKL)	
endif()

MACRO(HEADER_DIRECTORIES SearchPath return_list)
    FILE(GLOB_RECURSE new_list "${SearchPath}/*.h")
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

if (UseHYPRE)
	target_link_libraries(CESDSOL PRIVATE ${HYPRELibPath})
	HEADER_DIRECTORIES(${HYPREIncludePath} HYPREIncludeList)
	target_include_directories(CESDSOL PRIVATE ${HYPREIncludeList})
	target_compile_definitions(CESDSOL PRIVATE UseHYPRE)
endif()

target_include_directories(CESDSOL PRIVATE ${SourcesPath})

target_compile_definitions(CESDSOL PRIVATE MathLibrary=${MathBackend} ParallelismBackend=${Parallelism} ExamplesPath="${CMAKE_CURRENT_LIST_DIR}/examples/")
if (DebugMode)
	target_compile_definitions(CESDSOL PRIVATE DebugMode)
endif()
