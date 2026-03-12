add_executable(run_tests test_main.cpp)
target_link_libraries(run_tests PRIVATE fin_flutter_core)
add_test(NAME BenchmarkTests COMMAND run_tests)
