cmake_minimum_required(VERSION 3.22)
project(temporal_community_search)

find_package(TBB REQUIRED tbb)
add_library(temporal_community_search STATIC
        temporal_community_search.cpp
)

set_target_properties(temporal_community_search PROPERTIES PUBLIC_HEADER temporal_community_search.h)

target_link_libraries(temporal_community_search PUBLIC
        adj_list
        TBB::tbb
)

target_include_directories(temporal_community_search PUBLIC
        ../adj_list
)