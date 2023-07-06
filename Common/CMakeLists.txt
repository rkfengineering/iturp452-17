file(GLOB "LIB_SOURCES" src/*.cpp)
file(GLOB "LIB_HEADERS" include/*.h)

add_library(CommonLibrary STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(CommonLibrary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

#add_subdirectory(tests)
