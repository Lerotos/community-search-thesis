#include <gtest/gtest.h>
#include "temporal_community_search.h"
#include <cinttypes>

TEST(bigtest, testinputlarge){

    std::string filePath = "../../google_test/UnitTestInputs/TestInputLarge";
    temporal_community_search temporalCommunitySearch(filePath, 3, 2, 3);

    temporalCommunitySearch.TGR();

    std::set<int64_t> S{};
    auto t1bb = std::chrono::high_resolution_clock::now();

    temporalCommunitySearch.branchBound(S);

    auto t2bb = std::chrono::high_resolution_clock::now();
    auto ms_int_bb = std::chrono::duration_cast<std::chrono::milliseconds>(t2bb - t1bb);
    std::cout << "branchbound has taken " << ms_int_bb.count() << "ms\n" << std::endl;

    //this was not further checked, mostly used as comparison
    ASSERT_EQ(temporalCommunitySearch.R.size(), 3);
    ASSERT_EQ(temporalCommunitySearch.R.contains(4485), true);
    ASSERT_EQ(temporalCommunitySearch.R.contains(18220), true);
    ASSERT_EQ(temporalCommunitySearch.R.contains(45293), true);
}

TEST(massivetest, stackoverflowdata) {
    std::string filePath = "../../google_test/UnitTestInputs/stackoverflowdata_p1_mod";
    temporal_community_search temporalCommunitySearch(filePath, 150, 2, 200);

    temporalCommunitySearch.TGR();

    std::set<int64_t> S{};
    auto t1bb = std::chrono::high_resolution_clock::now();

    temporalCommunitySearch.branchBound(S);

    auto t2bb = std::chrono::high_resolution_clock::now();
    auto ms_int_bb = std::chrono::duration_cast<std::chrono::milliseconds>(t2bb - t1bb);
    std::cout << "branchbound has taken " << ms_int_bb.count() << "ms\n" << std::endl;

    printf("R: %" PRId64 "\n", temporalCommunitySearch.R.size());
}

TEST(massivetest_update, stackoverflowdata) {
    auto t1start = std::chrono::high_resolution_clock::now();

    std::string filePath = "../../google_test/UnitTestInputs/stackoverflowdata_p1";
    temporal_community_search temporalCommunitySearch(filePath, 150, 2, 200);
    printf("nr node aftewr inserting everything: %" PRId64 "\n", temporalCommunitySearch.adjList.vertices.size());

    auto t2start = std::chrono::high_resolution_clock::now();
    auto ms_int_start = std::chrono::duration_cast<std::chrono::milliseconds>(t2start - t1start);
    std::cout << "startup has taken " << ms_int_start.count() << "ms\n" << std::endl;



    auto t1tgr = std::chrono::high_resolution_clock::now();

    temporalCommunitySearch.TGR();
    printf("test %" PRId64 "\n", temporalCommunitySearch.adjList.neighbourMap.size());

    auto t2tgr = std::chrono::high_resolution_clock::now();
    auto ms_int_tgr = std::chrono::duration_cast<std::chrono::milliseconds>(t2tgr - t1tgr);
    std::cout << "tgr has taken " << ms_int_tgr.count() << "ms\n" << std::endl;


    auto t1u = std::chrono::high_resolution_clock::now();

    std::string filePath_update = "../../google_test/UnitTestInputs/stackoverflowdata_p2";
    temporalCommunitySearch.dynamic_update(filePath_update);
    printf("test %" PRId64 "\n", temporalCommunitySearch.adjList.neighbourMap.size());

    auto t2u = std::chrono::high_resolution_clock::now();
    auto ms_int_u = std::chrono::duration_cast<std::chrono::milliseconds>(t2u - t1u);
    std::cout << "update has taken " << ms_int_u.count() << "ms\n" << std::endl;

    printf("nr nodes after update: %" PRId64 "\n", temporalCommunitySearch.adjList.vertices.size());

    std::set<int64_t> S{};
    auto t1bb = std::chrono::high_resolution_clock::now();

    temporalCommunitySearch.branchBound(S);

    auto t2bb = std::chrono::high_resolution_clock::now();
    auto ms_int_bb = std::chrono::duration_cast<std::chrono::milliseconds>(t2bb - t1bb);
    std::cout << "branchbound has taken " << ms_int_bb.count() << "ms\n" << std::endl;
    printf("nr nodes after branchbound: %" PRId64 "\n", temporalCommunitySearch.adjList.vertices.size());


    printf("R: %" PRId64 "\n", temporalCommunitySearch.R.size());
}

TEST(massive_validationt_test, stackoverflowdata) {
    std::string filePath = "../../google_test/UnitTestInputs/paperinput";
    std::string filePath_update = "../../google_test/UnitTestInputs/stackoverflowdata_p1";
    std::string filePath_val = "../../google_test/UnitTestInputs/stackoverflowdata_p1_val";

    temporal_community_search temporalCommunitySearch(filePath, 150, 2, 200);

    temporalCommunitySearch.TGR();
    temporalCommunitySearch.reset();

    temporalCommunitySearch.dynamic_update(filePath_update);

    temporal_community_search validation(filePath_val, 150, 2, 200);

    validation.TGR();

    std::vector<int64_t> failureV;
    std::string DString = "Dv";
    std::string MIString = "MIv";

    auto func = [&](std::string& s){
        if (!failureV.empty()){
            int64_t sum = 0;
            for (auto n : failureV) {
                sum += n;
            }

            float avgOffset = static_cast<float>(sum) / static_cast<float>(failureV.size());

            std::cout << "Number of incorrect " << s << " values: " << failureV.size() << std::endl;
            std::cout << "Total offset: " << sum << std::endl;
            std::cout << "Average offset: " << avgOffset << std::endl;
        }

        failureV.clear();
    };

    for (auto& [n, rest]:temporalCommunitySearch.MICache) {
        auto MIv = temporalCommunitySearch.MICache[n];
        auto Dv = temporalCommunitySearch.DCache[n];
        auto MIvVal = validation.MICache[n];
        auto DvVal = validation.DCache[n];

        ASSERT_EQ(MIv.size(), Dv.size());
        ASSERT_EQ(MIv.size(), MIvVal.size());

        for (auto [start, end]:MIvVal) {
            ASSERT_EQ(end, MIv[start]) << failureV.emplace_back(std::abs(end - MIv[start]));
        }

        if (!failureV.empty()) std::cout << "Failure at node: " << n << std::endl;

        func(MIString);

        for (int i = 0; i < DvVal.size(); ++i) {
            ASSERT_EQ(DvVal[i], Dv[i]) << failureV.emplace_back(std::abs(DvVal[i] - Dv[i]));
        }

        if (!failureV.empty()) std::cout << "Failure at node: " << n << std::endl;

        func(DString);

        ASSERT_EQ(temporalCommunitySearch.degreeCache[n], validation.degreeCache[n])
                            << "degree offset at node : " << n << " is: "
                            << std::abs(temporalCommunitySearch.degreeCache[n] - validation.degreeCache[n])
                            << std::endl;
    }
}

TEST(tiny_insert_test, stackoverflowdata) {
    std::string filePath_update = "../../google_test/UnitTestInputs/paperinput";
    std::string filePath = "../../google_test/UnitTestInputs/stackoverflowdata_p1";

    temporal_community_search temporalCommunitySearch(filePath, 150, 2, 200);

    temporalCommunitySearch.TGR();

    temporalCommunitySearch.reset();


    temporalCommunitySearch.dynamic_update(filePath_update);

    temporal_community_search validation(filePath, 150, 2, 200);
    validation.TGR();
    validation.reset();

    validation.adjList.addFromFile(filePath_update);
    validation.TGR();

    std::vector<int64_t> failureV;
    int64_t failureD = 0;
    int64_t failureDval = 0;

    float Dsum = 0;
    int64_t nrOfDFailures = 0;

    std::vector<float> degreeFailures;
    int64_t cumDegree = 0;
    int64_t totalDegree = 0;

    auto funcMI = [&](){
        if (!failureV.empty()){
            int64_t sum = 0;
            for (auto n : failureV) {
                sum += n;
            }

            float avgOffset = static_cast<float>(sum) / static_cast<float>(failureV.size());

            std::cout << "Number of incorrect MIv values: " << failureV.size() << std::endl;
            std::cout << "Total offset: " << sum << std::endl;
            std::cout << "Average offset: " << avgOffset << std::endl;
        }

        failureV.clear();
    };

    auto funcD = [&](){
        if (failureD != 0){

            float avgFailRate = static_cast<float>(failureD) / static_cast<float>(failureDval);

            std::cout << "Nr of DCache failures: " << failureD << std::endl;
            std::cout << "Average fail rate of DCache: " << avgFailRate << std::endl;

            Dsum += avgFailRate;
            nrOfDFailures++;

        }

        failureD = 0;
        failureDval = 0;
    };

    for (auto& [n, rest]:temporalCommunitySearch.MICache) {
        auto MIv = temporalCommunitySearch.MICache[n];
        auto Dv = temporalCommunitySearch.DCache[n];
        auto MIvVal = validation.MICache[n];
        auto DvVal = validation.DCache[n];

        ASSERT_EQ(MIv.size(), Dv.size());
        ASSERT_EQ(MIv.size(), MIvVal.size());

        for (auto [start, end]:MIvVal) {
            EXPECT_EQ(end, MIv[start]) << failureV.emplace_back(std::abs(end - MIv[start]));
        }

        funcMI();

        for (int i = 0; i < DvVal.size(); ++i) {
            EXPECT_EQ(DvVal[i], Dv[i]) << failureD++;
            failureDval++;
        }

        funcD();


        EXPECT_EQ(temporalCommunitySearch.degreeCache[n], validation.degreeCache[n])
                            << degreeFailures.emplace_back(static_cast<float>(std::abs(temporalCommunitySearch.degreeCache[n]
                            - validation.degreeCache[n])));

        if (temporalCommunitySearch.degreeCache[n] != validation.degreeCache[n])
            totalDegree += validation.degreeCache[n];

        cumDegree += validation.degreeCache[n];
    }

    if (Dsum != 0){

        float totalFailRate = static_cast<float>(nrOfDFailures) / static_cast<float>(validation.DCache.size());
        float totalAvgFailRate = static_cast<float>(Dsum) / static_cast<float>(nrOfDFailures);

        std::cout << "Total DCache fail rate: " << totalFailRate << std::endl;
        std::cout << "Total DCache average fail rate : " << totalAvgFailRate << std::endl;


    }

    if (!degreeFailures.empty()){

        float sum = 0;
        for (auto degree:degreeFailures) {
            sum += degree;
        }

        float totalFailRate = static_cast<float>(degreeFailures.size()) / static_cast<float>(validation.degreeCache.size());
        float totalAvgFailRate = static_cast<float>(sum) / static_cast<float>(totalDegree);

        std::cout << "Total degree fail rate: " << totalFailRate << std::endl;
        std::cout << "Total degree average fail rate : " << totalAvgFailRate << std::endl;
    }
}