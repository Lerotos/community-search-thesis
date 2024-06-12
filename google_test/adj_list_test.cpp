#include <gtest/gtest.h>
#include "temporal_community_search.h"

TEST(computeComponents, testinput){
    std::string filePath = "../../google_test/UnitTestInputs/TestInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    std::unordered_set<int64_t> set{2, 3, 5, 6, 12, 15, 55, 66, 2364};

    auto components = temporalCommunitySearch.adjList.computeComponents(set);
    ASSERT_EQ(components.size(), 2);
    auto cLast = components.back();
    EXPECT_EQ(cLast.size(), 2);
    components.pop_back();
    auto cFirst = components.back();
    EXPECT_EQ(cFirst.size(), 7);

}

TEST(neighbourMap, testinput){
    std::string filePath = "../../google_test/UnitTestInputs/TestInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    auto neighbours = temporalCommunitySearch.adjList.neighbourMap;;

    ASSERT_EQ(neighbours.size(), 9);

    std::vector<int64_t> nodes{2, 3, 5, 6, 12, 15, 55, 66, 2364};
    std::vector<int64_t> neighbourSizes{3, 1, 1, 2, 1, 2, 1, 1, 2};


    for (int i = 0; i < nodes.size(); ++i) {
        EXPECT_EQ(neighbours[nodes[i]].size(), neighbourSizes[i]);
    }

    EXPECT_EQ(neighbours[2][6].size(), 5);
}

