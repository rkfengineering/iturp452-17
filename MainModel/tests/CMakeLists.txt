file(GLOB "TEST_SOURCES" *.cpp no_clutter_profile_tests/*.cpp)
file(GLOB "TEST_HEADERS" *.h no_clutter_profile_tests/*.h)

add_executable(
    ITUR_P452_test
    ${TEST_SOURCES}
    ${TEST_HEADERS}
)
target_link_libraries(
    ITUR_P452_test
    GTest::gtest_main
    CommonLibrary
    MainModel
    ClutterModel
)
include(GoogleTest)
gtest_discover_tests(ITUR_P452_test)