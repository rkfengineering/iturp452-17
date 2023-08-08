file(GLOB "TEST_SOURCES" *.cpp clutter_profile_tests/*.cpp)
file(GLOB "TEST_HEADERS" *.h clutter_profile_tests/*.h)

add_executable(
    ClutterModel_P452_17_test
    ${TEST_SOURCES}
    ${TEST_HEADERS}
)
target_link_libraries(
    ClutterModel_P452_17_test
    GTest::gtest_main
    CommonLibrary
    ITUR_P452_Model
    ClutterModel_P452_17
)

include(GoogleTest)
gtest_discover_tests(ClutterModel_P452_17_test)