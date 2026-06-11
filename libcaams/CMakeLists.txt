
file(GLOB MY_LIB_SOURCES CONFIGURE_DEPENDS "src/*.cpp" "src/*.c")
file(GLOB MY_LIB_HEADERS CONFIGURE_DEPENDS "include/*.h" "include/*.hpp")

# 3. Create the library target
add_library(caams STATIC ${MY_LIB_SOURCES} ${MY_LIB_HEADERS})

# 4. Link the include directory so files can find the headers
target_include_directories(caams PUBLIC include)
target_link_libraries(caams Eigen3::Eigen)
