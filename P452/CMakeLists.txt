file(GLOB "P452_SOURCES" src/*.cpp)
file(GLOB "P452_HEADERS" include/*.h)

add_library(P452Lib STATIC ${P452_SOURCES} ${P452_HEADERS})

target_include_directories(P452Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(P452Lib PUBLIC MainModel GasModel CommonLibrary ClutterModel)

add_subdirectory(tests)
