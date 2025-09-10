cmake_minimum_required(VERSION 3.13)

project(secsidh C ASM)

# Prepare variable to collect enabled parametersets
set(SECSIDH_PARAMETER_SETS_COMMON "")
set(SECSIDH_PARAMETER_SETS_CSIDH "")
set(SECSIDH_PARAMETER_SETS_CTIDH "")

# set(CMAKE_VERBOSE_MAKEFILE ON)

# define the options for this project
option(SECSIDH_GLOBAL_COUNTERS "Use global counters" ON)
option(SECSIDH_USE_OPENSSL_RAND "Use OpenSSL RNG" OFF)


option(SECSIDH_USE_GMP "Use the GMP implementation" ON)
option(SECSIDH_USE_GMP_SEC_FUNCTIONS "Use GMP low-level functions for cryptography" ON)
option(SECSIDH_USE_KARATSUBA "Use the Karatsuba multiplication" ON)
option(SECSIDH_USE_MONTGOMERY "Use Montgomery reduction" ON)
option(SECSIDH_USE_AVX2 "Use AVX2" OFF)

option(SECSIDH_USE_RANDOMNESS "Use the randomness-based implementation" OFF)
set(SECSIDH_CSIDH_VARIANT "df")

option(SECSIDH_CTIDH_TUNE_MULTS "Use version tuned for multiplications " ON)
option(SECSIDH_CTIDH_TUNE_CYCLES "Use version tuned for clock cycles " OFF)

option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)

# option(SECSIDH_USE_DUMMY "Use dummy isogenies" OFF)

# convenience macro define all parameter sets
macro(SECSIDH_DEFINE_PARAMETERS bits default)
    option(SECSIDH_ENABLE_CTIDH${bits} "Enable CTIDH${bits} parameterset." ${default})
    if(SECSIDH_ENABLE_CTIDH${bits})
        set(SECSIDH_PARAMETER_SETS_CTIDH ${SECSIDH_PARAMETER_SETS_CTIDH} ${bits})
    endif()
endmacro()
SECSIDH_DEFINE_PARAMETERS(2047m1l226 OFF)
SECSIDH_DEFINE_PARAMETERS(2047m4l205 ON)
SECSIDH_DEFINE_PARAMETERS(2047m6l194 ON)
SECSIDH_DEFINE_PARAMETERS(2047m7l188 OFF)


# handle counter setting
if(SECSIDH_GLOBAL_COUNTERS)
    add_definitions(-DGLOBAL_COUNTERS=1)
endif()

# if(SECSIDH_USE_DUMMY)
#     add_definitions(-DUSE_DUMMY)
# endif()

# Set project code options
# set(CMAKE_C_STANDARD 11)
# set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Pull in compiler flags
include(.CMake/compiler_options.cmake)

# Set up headers
configure_file(src/secsidhconfig.h.cmake
               ${PROJECT_BINARY_DIR}/include/secsidh/secsidhconfig.h)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/secsidh.h ${CMAKE_BINARY_DIR}/include/secsidh/secsidh.h)
include_directories(${PROJECT_BINARY_DIR}/include)

# Set subdirectory for sources
add_subdirectory(src)

# Set up library
foreach(bits IN ITEMS ${SECSIDH_PARAMETER_SETS_CTIDH})
    if(SECSIDH_ENABLE_CTIDH${bits})
        set(SOURCES ${SOURCES}
            $<TARGET_OBJECTS:secsidh-ctidh${bits}>
            $<TARGET_OBJECTS:csidh-ctidh${bits}>
            $<TARGET_OBJECTS:ctidh-${bits}>
            $<TARGET_OBJECTS:uint-${bits}>
            $<TARGET_OBJECTS:primes-ctidh${bits}>
            $<TARGET_OBJECTS:fp-${bits}>
            $<TARGET_OBJECTS:poly-ctidh${bits}>
            $<TARGET_OBJECTS:mont-ctidh${bits}>
            $<TARGET_OBJECTS:elligator-ctidh${bits}>
            $<TARGET_OBJECTS:skgen-ctidh${bits}>
            $<TARGET_OBJECTS:validate-ctidh${bits}>
            $<TARGET_OBJECTS:csidh-ctidh${bits}>
            $<TARGET_OBJECTS:steps-ctidh${bits}>
            $<TARGET_OBJECTS:random-ctidh${bits}>
            $<TARGET_OBJECTS:steps_t-ctidh${bits}>
            $<TARGET_OBJECTS:int32_sort>
            $<TARGET_OBJECTS:fips202>
            $<TARGET_OBJECTS:keccakf1600>
        )
        if(SECSIDH_USE_KARATSUBA)
            set(SOURCES ${SOURCES} $<TARGET_OBJECTS:karatsuba_fp_sqr>)
        endif()
    endif()
endforeach()


add_library(secsidh STATIC ${SOURCES})
if(SECSIDH_USE_GMP)
    target_link_libraries(secsidh PUBLIC gmp rng)
else()
   target_link_libraries(secsidh PUBLIC rng)
endif()



# Main executables (XXX: rename to examples)
add_subdirectory(main)

set_target_properties(secsidh
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    VERSION "0.1"
    SOVERSION 0
)

install(TARGETS secsidh
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
)

install(FILES
        ${PROJECT_SOURCE_DIR}/src/secsidh.h
        ${PROJECT_BINARY_DIR}/include/secsidh/secsidhconfig.h
        DESTINATION include/secsidh
)

# Add benchmark targets
option(SECSIDH_BENCHMARK_RUNS "Number of runs for each benchmark" 100)

# Ensure BENCHMARK_RUNS is a number
if(NOT SECSIDH_BENCHMARK_RUNS OR SECSIDH_BENCHMARK_RUNS STREQUAL "OFF")
    set(SECSIDH_BENCHMARK_RUNS 100)
    message(STATUS "Setting SECSIDH_BENCHMARK_RUNS to default value: 100")
endif()

# Create benchmark target for CTIDH
foreach(bits IN ITEMS ${SECSIDH_PARAMETER_SETS_CTIDH})
    if(SECSIDH_ENABLE_CTIDH${bits})
        add_custom_target(benchmark-ctidh-${bits}
            COMMAND ${CMAKE_COMMAND} -E echo "Running CTIDH-${bits} benchmark with ${SECSIDH_BENCHMARK_RUNS} iterations..."
            COMMAND ${CMAKE_BINARY_DIR}/main/ctidh-${bits}.main -bact ${SECSIDH_BENCHMARK_RUNS} > ${CMAKE_BINARY_DIR}/benchmark-ctidh-${bits}.log
            COMMAND ${CMAKE_COMMAND} -E echo "Analyzing benchmark results for CTIDH-${bits}..."
            COMMAND python3 ${CMAKE_SOURCE_DIR}/analyze_bench.py < ${CMAKE_BINARY_DIR}/benchmark-ctidh-${bits}.log > ${CMAKE_BINARY_DIR}/benchmark-ctidh-${bits}-analysis.log
            COMMAND ${CMAKE_COMMAND} -E echo "=== BENCHMARK RESULTS FOR CTIDH-${bits} ==="
            COMMAND cat ${CMAKE_BINARY_DIR}/benchmark-ctidh-${bits}-analysis.log
            COMMAND ${CMAKE_COMMAND} -E echo "=========================================="
            DEPENDS ctidh-${bits}.main
            COMMENT "Benchmarking CTIDH-${bits}"
            VERBATIM
        )
    endif()
endforeach()

# Add a benchmark summary target that does not depend on benchmark-all
add_custom_target(benchmark-summary
    COMMAND ${CMAKE_COMMAND} -E echo "=== BENCHMARK SUMMARY ==="
    COMMENT "Showing benchmark summary"
    VERBATIM
)

# Add commands to print benchmark summaries for each parameter set
foreach(bits IN ITEMS ${SECSIDH_PARAMETER_SETS_CTIDH})
    if(SECSIDH_ENABLE_CTIDH${bits})
        add_custom_command(
            TARGET benchmark-summary
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "Results for CTIDH-${bits}:"
            COMMAND grep -A 20 "=== ACTION ===" ${CMAKE_BINARY_DIR}/benchmark-ctidh-${bits}-analysis.log || echo "No action results found"
            COMMAND ${CMAKE_COMMAND} -E echo "------------------------------------------"
            VERBATIM
        )
    endif()
endforeach()

# Master benchmark target that just runs the benchmarks (no summary)
add_custom_target(benchmark-all
    COMMENT "Running all benchmarks"
)

# Add dependencies to the benchmark-all target
foreach(bits IN ITEMS ${SECSIDH_PARAMETER_SETS_CTIDH})
    if(SECSIDH_ENABLE_CTIDH${bits})
        add_dependencies(benchmark-all benchmark-ctidh-${bits})
    endif()
endforeach()

# Combined target that runs benchmarks and then shows the summary
add_custom_target(benchmark
    DEPENDS benchmark-all benchmark-summary
    COMMENT "Running all benchmarks and showing summary"
)
