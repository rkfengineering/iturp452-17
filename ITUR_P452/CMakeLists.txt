file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(ITUR_P452_Model STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(ITUR_P452_Model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ITUR_P452_Model LINK_PUBLIC ClutterModel_P452_17 CommonLibrary GasModel GTest::gtest_main)

add_subdirectory(tests)
