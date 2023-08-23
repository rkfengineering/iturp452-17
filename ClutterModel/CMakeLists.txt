file(GLOB "CLUTTER_SOURCES" src/*.cpp)
file(GLOB "CLUTTER_HEADERS" include/*.h)

add_library(ClutterModel STATIC ${CLUTTER_SOURCES} ${CLUTTER_HEADERS})

target_include_directories(ClutterModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ClutterModel LINK_PUBLIC MainModel CommonLibrary GTest::gtest_main)

add_subdirectory(tests)
