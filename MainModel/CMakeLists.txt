file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(MainModel STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(MainModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MainModel LINK_PUBLIC ClutterModel CommonLibrary GasModel GTest::gtest_main)

add_subdirectory(tests)
