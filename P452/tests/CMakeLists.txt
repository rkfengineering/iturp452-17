file(GLOB "TEST_SOURCES" *.cpp)
file(GLOB "TEST_HEADERS" *.h)

add_executable(
    P452_wrapper_test
    ${TEST_SOURCES}
    ${TEST_HEADERS}
)
target_link_libraries(
    P452_wrapper_test
    GTest::gtest_main
    P452Lib
)
include(GoogleTest)
gtest_discover_tests(P452_wrapper_test)