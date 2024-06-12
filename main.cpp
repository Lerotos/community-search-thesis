#include "adj_list.h"
#include "temporal_community_search.h"
#include <cinttypes>

int main() {
    std::string filePath = "../testinput";

    //temporal_community_search temporalCommunitySearch(filePath, 3, 2);


    /*
     * for testing purposes
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
     */


    return 0;
}
