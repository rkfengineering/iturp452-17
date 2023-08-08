file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(ClutterModel_P452_17 STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(ClutterModel_P452_17 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ClutterModel_P452_17 LINK_PUBLIC ITUR_P452_Model CommonLibrary GTest::gtest_main)

add_subdirectory(tests)
