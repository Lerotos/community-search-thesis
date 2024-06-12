#ifndef TEMPUS_ADJ_LIST_H
#define TEMPUS_ADJ_LIST_H

#include <set>
#include <map>
#include <cstdint>
#include <unordered_set>
#include "dynamic_graph/dynamic_connectivity.hpp"

class AdjList{
private:

    //graph size has to be predetermined as DynamicConnectivity does not yet support graph resizing.
    const int64_t k = 300000;

    /**
     * Computes the extracted information from addFromFile and inserts or removes them from the graph
     * and the edge set accordingly after checking if they already exist.
     * Deletes edges after insertions. Meaning inserting -> deleting -> reinserting the same node
     * results in the node not being added.
     *
     * @param adds
     * @param dels
     */
    void computeBatch(std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>> &addMap,
                      std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>> &delMap);

public:
    DynamicConnectivity dynamicConnectivity;

    //keeps track of all vertices currently in the graph
    std::unordered_set<int64_t> vertices;

    //                 vertex                      neighbour         timestamps of the edges
    std::unordered_map<int64_t, std::unordered_map<int64_t, std::set<int64_t>>> neighbourMap;

    AdjList();

    /**
     * Used to compute the connected components of the given graph.
     *
     * @param uniqueVerticesSet set of unique nodes
     * @return vector of component sets
     */
    std::vector<std::set<int64_t>>
    computeComponents(std::unordered_set<int64_t> &uniqueVerticesSet) const;

    /**
     * Reads the given file and extracts information over nodes that are to be added or removed from the graph.
     *
     * @param path to the file
     */
    std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>>
    addFromFile(const std::string &path);

};

#endif //TEMPUS_ADJ_LIST_H
