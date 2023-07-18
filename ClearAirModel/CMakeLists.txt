file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(ClearAirModel STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(ClearAirModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ClearAirModel LINK_PUBLIC CommonLibrary GasModel)

add_subdirectory(tests)
