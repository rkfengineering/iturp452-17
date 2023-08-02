file(GLOB "TEST_SOURCES" *.cpp clutter_profile_tests/*.cpp)
file(GLOB "TEST_HEADERS" *.h clutter_profile_tests/*.h)

add_executable(
    ClearAirModel_test
    ${TEST_SOURCES}
    ${TEST_HEADERS}
)
target_link_libraries(
    ClearAirModel_test
    GTest::gtest_main
    CommonLibrary
    ClearAirModel
)

include(GoogleTest)
gtest_discover_tests(ClearAirModel_test)