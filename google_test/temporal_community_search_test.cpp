#include <gtest/gtest.h>
#include "temporal_community_search.h"
#include <cinttypes>


TEST(genVertexVector, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    auto Nu = temporalCommunitySearch.adjList.neighbourMap[1];
    auto vertexVector = temporalCommunitySearch.genVertexVector(Nu);

    std::unordered_set<int64_t> set{2, 3, 4};
    ASSERT_EQ(vertexVector.size(), 3);
    for (int64_t v:vertexVector) {
        ASSERT_EQ(set.contains(v), true);
    }
}

TEST(meta_interval_decomposition, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.meta_interval_decomposition(4);
    auto MIu = temporalCommunitySearch.MICache[4];

    ASSERT_EQ(MIu.size(), 3);
    EXPECT_EQ(MIu[6], 8);
    EXPECT_EQ(MIu[8], 9);
    EXPECT_EQ(MIu[9], 11);

    std::vector testVector{1, 2, 1};
    auto Du = temporalCommunitySearch.DCache[4];

    ASSERT_EQ(Du.size(), 3);
    for (int i = 0; i < Du.size(); ++i) {
        EXPECT_EQ(Du[i], testVector[i]);
    }
}

TEST(meta_interval_decomposition_repeated_neighbours, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.meta_interval_decomposition(1);
    auto MIu = temporalCommunitySearch.MICache[1];
    ASSERT_EQ(MIu.size(), 7);
    EXPECT_EQ(MIu[1], 2);
    EXPECT_EQ(MIu[2], 3);
    EXPECT_EQ(MIu[3], 5);
    EXPECT_EQ(MIu[5], 6);
    EXPECT_EQ(MIu[6], 7);
    EXPECT_EQ(MIu[7], 9);
    EXPECT_EQ(MIu[9], 10);

    std::vector testVector{1, 2, 2, 1, 1, 2, 1};

    auto Du = temporalCommunitySearch.DCache[1];
    ASSERT_EQ(Du.size(), 7);
    for (int i = 0; i < Du.size(); ++i) {
        EXPECT_EQ(Du[i], testVector[i]);
    }
}

TEST(deleteSingleNode, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), true);

    ASSERT_EQ(temporalCommunitySearch.deletedNodes.find(1) != temporalCommunitySearch.deletedNodes.end(), false);

    std::vector<UndirectedEdge> deletionCache;
    temporalCommunitySearch.deleteSingleNode(1, deletionCache);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), false);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), false);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), false);

    ASSERT_EQ(temporalCommunitySearch.deletedNodes.find(1) != temporalCommunitySearch.deletedNodes.end(), true);

}

TEST(deleteNodes, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), true);

    ASSERT_EQ(temporalCommunitySearch.deletedNodes.find(1) != temporalCommunitySearch.deletedNodes.end(), false);

    std::vector<int64_t> A{1};
    std::vector<UndirectedEdge> deletionCache = temporalCommunitySearch.deleteNodes(A);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), false);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), false);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), false);

    ASSERT_EQ(temporalCommunitySearch.deletedNodes.find(1) != temporalCommunitySearch.deletedNodes.end(), true);

}

TEST(reAddNode, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    std::vector<UndirectedEdge> deletionCache;
    temporalCommunitySearch.deleteSingleNode(1, deletionCache);

    temporalCommunitySearch.reAddNodes(deletionCache);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), true);

    ASSERT_EQ(temporalCommunitySearch.deletedNodes.find(1) != temporalCommunitySearch.deletedNodes.end(), false);
}

TEST(computeDegreePersistence, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.meta_interval_decomposition(1);
    temporalCommunitySearch.meta_interval_decomposition(2);
    temporalCommunitySearch.meta_interval_decomposition(3);
    temporalCommunitySearch.meta_interval_decomposition(4);

    EXPECT_EQ(temporalCommunitySearch.computeDegreePersistence(1), 8);
    EXPECT_EQ(temporalCommunitySearch.computeDegreePersistence(2), 5);
    EXPECT_EQ(temporalCommunitySearch.computeDegreePersistence(3), 6);
    EXPECT_EQ(temporalCommunitySearch.computeDegreePersistence(4), 4);
}

TEST(TGR, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), true);

    temporalCommunitySearch.TGR();

    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 2), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 3), true);
    ASSERT_EQ(temporalCommunitySearch.adjList.dynamicConnectivity.IsConnected(1, 4), true);

    temporal_community_search temporalCommunitySearch2(filePath, 3, 2, 5);
    temporalCommunitySearch2.TGR();

    ASSERT_EQ(temporalCommunitySearch2.adjList.dynamicConnectivity.IsConnected(1, 2), false);
    ASSERT_EQ(temporalCommunitySearch2.adjList.dynamicConnectivity.IsConnected(1, 3), false);
    ASSERT_EQ(temporalCommunitySearch2.adjList.dynamicConnectivity.IsConnected(1, 4), false);


}

TEST(removeNode, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 5);
    std::set<int64_t> set{2};

    temporalCommunitySearch.meta_interval_decomposition(1);
    temporalCommunitySearch.meta_interval_decomposition(2);
    temporalCommunitySearch.meta_interval_decomposition(3);
    temporalCommunitySearch.meta_interval_decomposition(4);
    temporalCommunitySearch.degreeCache[4] = 0;
    temporalCommunitySearch.degreeCache[1] = 8;
    temporalCommunitySearch.degreeCache[2] = 5;
    temporalCommunitySearch.degreeCache[3] = 6;

    auto data = temporalCommunitySearch.removeNode(4, set);

    bool flag = data.has_value();
    auto payload = data.value();
    ASSERT_EQ(flag, true);
    ASSERT_EQ(payload.size(), 1);

    auto data2 = temporalCommunitySearch.removeNode(3, set);
    bool flag2 = data2.has_value();
    ASSERT_EQ(flag2, false);
}

TEST(addNode, paperinput) {
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 5);
    std::set<int64_t> set{};

    temporalCommunitySearch.meta_interval_decomposition(1);
    temporalCommunitySearch.meta_interval_decomposition(2);
    temporalCommunitySearch.meta_interval_decomposition(3);
    temporalCommunitySearch.meta_interval_decomposition(4);
    temporalCommunitySearch.degreeCache[4] = 0;
    temporalCommunitySearch.degreeCache[1] = 8;
    temporalCommunitySearch.degreeCache[2] = 5;
    temporalCommunitySearch.degreeCache[3] = 6;

    auto data = temporalCommunitySearch.removeNode(4, set);

    ASSERT_EQ(temporalCommunitySearch.degreeCache[3], 4);

    temporalCommunitySearch.addNode(4);

    ASSERT_EQ(temporalCommunitySearch.degreeCache[3], 6);
    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);
}

TEST(meta_interval_intersection, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.meta_interval_decomposition(1);
    temporalCommunitySearch.meta_interval_decomposition(2);
    temporalCommunitySearch.meta_interval_decomposition(3);
    temporalCommunitySearch.meta_interval_decomposition(4);

    std::set<int64_t> set0{};
    std::set<int64_t> set1{1};
    std::set<int64_t> set2{1,2};
    std::set<int64_t> set3{1,2,3};
    std::set<int64_t> set4{1,2,3,4};

    auto result1 = temporalCommunitySearch.meta_interval_intersection(set1);
    ASSERT_EQ(result1, 8);

    auto result2 = temporalCommunitySearch.meta_interval_intersection(set2);

    ASSERT_EQ(result2, 4);

    auto result3 = temporalCommunitySearch.meta_interval_intersection(set3);

    ASSERT_EQ(result3, 4);

    auto result4 = temporalCommunitySearch.meta_interval_intersection(set4);

    ASSERT_EQ(result4, 0);

    auto result0 = temporalCommunitySearch.meta_interval_intersection(set0);

    ASSERT_EQ(result0, std::numeric_limits<int64_t>::max());
}

TEST(branchbound, paperinput) {
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    std::set<int64_t> S{};

    temporalCommunitySearch.branchBound(S);

    ASSERT_EQ(temporalCommunitySearch.R.size(), 3);
    ASSERT_FALSE(temporalCommunitySearch.R.find(1) == temporalCommunitySearch.R.end());
    ASSERT_FALSE(temporalCommunitySearch.R.find(3) == temporalCommunitySearch.R.end());
    ASSERT_FALSE(temporalCommunitySearch.R.find(2) == temporalCommunitySearch.R.end());
}




TEST(meta_interval_decomposition_update_1, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    auto MIC = temporalCommunitySearch.MICache[1];
    auto DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v{1, 2, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 7);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v.size());

    for (int i = 0; i < v.size(); ++i) {
        ASSERT_EQ(DC[i], v[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);

    std::set<int64_t> set {};
    temporalCommunitySearch.meta_interval_decomposition_update(1, 0, 20, set);

    MIC = temporalCommunitySearch.MICache[1];
    DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v2{1, 2, 3, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 8);

    ASSERT_EQ(MIC[0], 1);
    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v2.size());


    for (int i = 0; i < v2.size(); ++i) {
        ASSERT_EQ(DC[i], v2[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 9);

}

TEST(meta_interval_decomposition_update_2, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    auto MIC = temporalCommunitySearch.MICache[1];
    auto DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v{1, 2, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 7);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v.size());

    for (int i = 0; i < v.size(); ++i) {
        ASSERT_EQ(DC[i], v[i]);
    }


    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);

    std::set<int64_t> set {};
    temporalCommunitySearch.meta_interval_decomposition_update(1, 5, 20, set);

    MIC = temporalCommunitySearch.MICache[1];
    DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v2{1, 2, 2, 2, 2, 3, 2, 1};

    ASSERT_EQ(MIC.size(), 8);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 8);
    ASSERT_EQ(MIC[8], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v2.size());


    for (int i = 0; i < v2.size(); ++i) {
        ASSERT_EQ(DC[i], v2[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 10);
}

TEST(meta_interval_decomposition_update_3, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    auto MIC = temporalCommunitySearch.MICache[1];
    auto DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v{1, 2, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 7);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v.size());

    for (int i = 0; i < v.size(); ++i) {
        ASSERT_EQ(DC[i], v[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);


    std::set<int64_t> set {};

    temporalCommunitySearch.meta_interval_decomposition_update(1, 8, 20, set);

    MIC = temporalCommunitySearch.MICache[1];
    DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v2{1, 2, 2, 1, 1, 2, 3, 2, 1};

    ASSERT_EQ(MIC.size(), 9);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 8);
    ASSERT_EQ(MIC[8], 9);
    ASSERT_EQ(MIC[9], 10);
    ASSERT_EQ(MIC[10], 11);

    ASSERT_EQ(DC.size(), v2.size());


    for (int i = 0; i < v2.size(); ++i) {
        ASSERT_EQ(DC[i], v2[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 9);
}

TEST(meta_interval_decomposition_update_4, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    auto MIC = temporalCommunitySearch.MICache[1];
    auto DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v{1, 2, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 7);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v.size());

    for (int i = 0; i < v.size(); ++i) {
        ASSERT_EQ(DC[i], v[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);

    std::set<int64_t> set {};

    temporalCommunitySearch.meta_interval_decomposition_update(1, 20, 20, set);

    MIC = temporalCommunitySearch.MICache[1];
    DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v2{1, 2, 2, 1, 1, 2, 1, 1};

    ASSERT_EQ(MIC.size(), 8);

    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);
    ASSERT_EQ(MIC[20], 23);

    ASSERT_EQ(DC.size(), v2.size());


    for (int i = 0; i < v2.size(); ++i) {
        ASSERT_EQ(DC[i], v2[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);

}

TEST(dynamic_update, paperinput){
    std::string filePath = "../../google_test/UnitTestInputs/PaperInput";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    std::string updatePath = "../../google_test/UnitTestInputs/updateInput";
    temporalCommunitySearch.dynamic_update(updatePath);

    auto MIC = temporalCommunitySearch.MICache[1];
    auto DC = temporalCommunitySearch.DCache[1];
    std::vector<int64_t> v2{0, 1, 2, 2, 1, 1, 2, 1};

    ASSERT_EQ(MIC.size(), 8);

    ASSERT_EQ(MIC[0], 1);
    ASSERT_EQ(MIC[1], 2);
    ASSERT_EQ(MIC[2], 3);
    ASSERT_EQ(MIC[3], 5);
    ASSERT_EQ(MIC[5], 6);
    ASSERT_EQ(MIC[6], 7);
    ASSERT_EQ(MIC[7], 9);
    ASSERT_EQ(MIC[9], 10);

    ASSERT_EQ(DC.size(), v2.size());

    for (int i = 0; i < v2.size(); ++i) {
        ASSERT_EQ(DC[i], v2[i]);
    }

    ASSERT_EQ(temporalCommunitySearch.degreeCache[1], 8);

}

TEST(insert_paperInput, input1) {

    std::string filePath = "../../google_test/UnitTestInputs/input1";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();
    temporalCommunitySearch.reset();

    std::string filePath_update = "../../google_test/UnitTestInputs/paperinput";
    temporalCommunitySearch.dynamic_update(filePath_update);

    std::map<int64_t, int64_t> verification_map_1 {{1, 2}, {2, 3}, {3, 5}, {5, 6}, {6, 7}, {7, 9}, {9, 10}, {6666, 6669}};
    std::map<int64_t, int64_t> verification_map_2 {{1, 3}, {3, 4}, {4, 6}, {6, 7}, {8888, 8891}};
    std::map<int64_t, int64_t> verification_map_3 {{2, 4}, {4, 5}, {5, 7}, {7, 8}, {8, 10}, {10, 11}, {666, 669}};
    std::map<int64_t, int64_t> verification_map_4 {{6, 8}, {8, 9}, {9, 11}, {66, 69}};
    std::vector<std::map<int64_t, int64_t>> mapVector {verification_map_1, verification_map_2, verification_map_3, verification_map_4};

    for (int i = 1; i < 5; ++i) {
        auto controlMap = mapVector[i - 1];
        auto MIv = temporalCommunitySearch.MICache[i];
        ASSERT_EQ(MIv.size(), controlMap.size());

        for (auto [start, end]:MIv) {
            ASSERT_EQ(end, controlMap[start]);
        }
    }

    std::vector<int64_t> verification_vector_1 {1, 2, 2, 1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_2 {1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_3 {1, 2, 1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_4 {1, 2, 1, 0};
    std::vector<std::vector<int64_t>> vectorVector {verification_vector_1, verification_vector_2, verification_vector_3, verification_vector_4};

    for (int i = 0; i < 4; ++i) {
        auto controlVector = vectorVector[i];
        auto Dv = temporalCommunitySearch.DCache[i + 1];
        ASSERT_EQ(Dv.size(), mapVector[i].size());
        ASSERT_EQ(Dv.size(), controlVector.size());

        for (int j = 0; j < Dv.size(); ++j) {
            ASSERT_EQ(Dv[j], controlVector[j]);
        }
    }

    auto degrees = temporalCommunitySearch.degreeCache;
    ASSERT_EQ(degrees[1], 8);
    ASSERT_EQ(degrees[2], 5);
    ASSERT_EQ(degrees[3], 6);
    ASSERT_EQ(degrees[4], 4);
}

TEST(insert_paperInput_reverse, input1) {

    std::string filePath = "../../google_test/UnitTestInputs/input1";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();
    temporalCommunitySearch.reset();

    std::string filePath_update = "../../google_test/UnitTestInputs/paperinputreverse";
    temporalCommunitySearch.dynamic_update(filePath_update);

    for (auto& xx:temporalCommunitySearch.MICache) {
        printf("MI_%" PRId64 " = {", xx.first);
        for (auto & [start, end] : xx.second) {
            printf("(%" PRId64 ", %" PRId64 "), ", start, end);
        }
        printf("}\n");
    }

    for (auto& xxx:temporalCommunitySearch.DCache) {
        printf("D_%" PRId64 " = {", xxx.first);
        for (auto x : xxx.second) {
            printf("%" PRId64 ", ", x);
        }
        printf("}\n");
    }
    printf("}\n");

    std::map<int64_t, int64_t> verification_map_1 {{1, 2}, {2, 3}, {3, 5}, {5, 6}, {6, 7}, {7, 9}, {9, 10}, {6666, 6669}};
    std::map<int64_t, int64_t> verification_map_2 {{1, 3}, {3, 4}, {4, 6}, {6, 7}, {8888, 8891}};
    std::map<int64_t, int64_t> verification_map_3 {{2, 4}, {4, 5}, {5, 7}, {7, 8}, {8, 10}, {10, 11}, {666, 669}};
    std::map<int64_t, int64_t> verification_map_4 {{6, 8}, {8, 9}, {9, 11}, {66, 69}};
    std::vector<std::map<int64_t, int64_t>> mapVector {verification_map_1, verification_map_2, verification_map_3, verification_map_4};

    for (int i = 1; i < 5; ++i) {
        auto controlMap = mapVector[i - 1];
        auto MIv = temporalCommunitySearch.MICache[i];
        ASSERT_EQ(MIv.size(), controlMap.size());

        for (auto [start, end]:MIv) {
            ASSERT_EQ(end, controlMap[start]);
        }
    }

    std::vector<int64_t> verification_vector_1 {1, 2, 2, 1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_2 {1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_3 {1, 2, 1, 1, 2, 1, 0};
    std::vector<int64_t> verification_vector_4 {1, 2, 1, 0};
    std::vector<std::vector<int64_t>> vectorVector {verification_vector_1, verification_vector_2, verification_vector_3, verification_vector_4};

    for (int i = 0; i < 4; ++i) {
        auto controlVector = vectorVector[i];
        auto Dv = temporalCommunitySearch.DCache[i + 1];
        ASSERT_EQ(Dv.size(), mapVector[i].size());
        ASSERT_EQ(Dv.size(), controlVector.size());

        for (int j = 0; j < Dv.size(); ++j) {
            ASSERT_EQ(Dv[j], controlVector[j]);
        }
    }

    auto degrees = temporalCommunitySearch.degreeCache;
    ASSERT_EQ(degrees[1], 8);
    ASSERT_EQ(degrees[2], 5);
    ASSERT_EQ(degrees[3], 6);
    ASSERT_EQ(degrees[4], 4);
}