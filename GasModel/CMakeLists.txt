file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(GasModel STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(GasModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(GasModel LINK_PUBLIC CommonLibrary)

add_subdirectory(tests)
