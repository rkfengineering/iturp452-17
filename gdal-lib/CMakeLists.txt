file(GLOB "GDAL_SOURCES" src/*.cpp)
file(GLOB "GDAL_HEADERS" include/gdal_lib/*.h)

add_library(GdalLib STATIC ${GDAL_SOURCES} ${GDAL_HEADERS})

target_include_directories(GdalLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(GdalLib PUBLIC GDAL::GDAL)

#add_subdirectory(tests)