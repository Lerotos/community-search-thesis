cmake_minimum_required(VERSION 3.22)
project(tempus)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(tempus main.cpp)

#add_subdirectory(include/libcuckoo)
#add_subdirectory(include/parlaylib/)
add_subdirectory(include/adj_list)
add_subdirectory(include/temporal_community_search)
add_subdirectory(include/dynamic-connectivity-hdt/src/dynamic_graph)
add_subdirectory(include/dynamic-connectivity-hdt/src/utilities)
add_subdirectory(google_test)

target_link_libraries(tempus PUBLIC
        temporal_community_search
)

target_include_directories(tempus PUBLIC
        include/temporal_community_search
)

