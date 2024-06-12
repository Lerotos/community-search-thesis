#ifndef TEMPUS_TEMPORAL_COMMUNITY_SEARCH_H
#define TEMPUS_TEMPORAL_COMMUNITY_SEARCH_H


#include <cstdint>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>
#include "dynamic_graph/dynamic_connectivity.hpp"
#include "adj_list.h"
#include <tbb/tbb.h>

class temporal_community_search {

public:
    int64_t size = 0;
    int64_t theta;
    int64_t k;
    int64_t tau;

    tbb::concurrent_unordered_map<int64_t, std::map<int64_t, int64_t>> MICache;
    tbb::concurrent_unordered_map<int64_t, std::vector<int64_t>> DCache;
    tbb::concurrent_unordered_map<int64_t, int64_t> degreeCache;

    std::set<int64_t> R;
    std::unordered_set<int64_t> deletedNodes;
    std::unordered_set<int64_t> preCutNodes;
    std::vector<UndirectedEdge> TGRDeletionCache;

    AdjList adjList;

    /**
     * Initialises the initial graph, theta and k
     *
     * @param path to the file containing the edges to be added
     * @param theta length of the time interval
     * @param k core requirement
     * @param tau persistence threshold
     */
    explicit temporal_community_search(std::string &path, int64_t theta, int64_t k, int64_t tau);

    /**
     * Algorithm 1 from the paper. Generates the meta interval array MI and the 𝜃-persistent degree array D.
     *
     * @param u given vertex
     */
    void meta_interval_decomposition(int64_t u);

    /**
     * Algorithm 2 from the paper. Preprocessing by pruning the graph of vertices that fall below the persistence
     * threshold.
     */
    void TGR();

    /**
     * Computes the degree persistence of the given node.
     *
     * @param u given node
     * @return degree persistence
     */
    int64_t computeDegreePersistence(int64_t u);


    /**
     * Algorithm 3 from the paper. Calculates the common meta intervals of the given nodes.
     *
     * @param components set of connected nodes
     * @return the total length + theta
     */
    int64_t meta_interval_intersection(std::set<int64_t> &components);

    /**
     * Algorithm 4 from the paper. changes the persistence values of the nodes neighbours (by using the
     * removeNodeHelper function) and adds nodes to the queue if they need to be removed. If any of these nodes are
     * part of @param S the function will return empty.
     *
     * @param u given node
     * @param S set of nodes must not be removed
     * @return either a set of nodes to be removed or an empty value if any of them are in S
     */
    std::optional<std::unordered_set<int64_t>> removeNode(int64_t u, std::set<int64_t> &S);

    /**
     * Algorithm 5 from the paper. Used to revert the changes made by the removeNode algorithm.
     *
     * @param u given node
     */
    void addNode(int64_t u);

    /**
     * Branch-Bound Algorithm from the paper. Used to find the largest persistent k-core.
     *
     * @param S set of nodes that must be included in the result set
     */
    void branchBound(std::set<int64_t> &S);

    /**
     * Adds edges from the @param deletionCache back into the graph
     *
     * @param deletionCache
     */
    void reAddNodes(std::vector<UndirectedEdge> &deletionCache);

    /**
     * Deletes the vertices in @param A from the graph and stores them in a vector which is then returned.
     * Uses deleteSingleNodes to do so.
     *
     * @param removeVector vertices that are to be deleted
     * @return Cache of deleted vertices
     */
    std::vector<UndirectedEdge> deleteNodes(std::vector<int64_t> &removeVector);


    /**
     * Deletes a single node from the graph.
     *
     * @param u node to be deleted
     * @param deletionCache storage for the deleted edges
     */
    void deleteSingleNode(int64_t u, std::vector<UndirectedEdge> &deletionCache);

    /**
     * Helper function for TGR and removeNode to avoid duplicate code.
     * Iterates through the given nodes neighbours and updates their persistence values if needed.
     * Calls the function @param f which are the differences between removeNode and TGR.
     *
     * @tparam F function parameter
     * @param u given node
     * @param f function
     */
    template<typename F>
    void removeNodeHelper(int64_t u, F &&f);

    /**
     * Helper function for branchBound to avoid duplicate code.
     * Responsible for the removeNode calls as well as recursion.
     *
     * @param Q Queue is will be iterated over
     * @param A vector to remember all nodes that were iterated over
     * @param S saved nodes
     */
    void branchBoundHelper(std::queue<int64_t> &Q, std::set<int64_t> &S);

    /**
     * Generates a vector of active neighbours.
     *
     * @param Nu Map containing the neighbours (and their timestamps)
     * @return the generated vector
     */
    template<typename T>
    std::vector<int64_t> genVertexVector(const T& Nu) const;

    /**
     * Resets the graph to the state before pruning was applied and resets size and R.
     */
    void reset();

    /**
     * Adjusts MICache, DCache and degreeCache after an insert without having to recalculate from scratch.
     *
     * @param u given node
     * @param time of the edge
     * @param neighbour of the given node of that edge
     * @param neighbourTimes set of times where the neighbour appeared before
     */
    void meta_interval_decomposition_update(int64_t u, int64_t time, int64_t neighbour, std::set<int64_t> &neighbourTimes);

    /**
     * Dynamically integrates the edges given through the filepath.
     *
     * @param path to the file containing the edges to be added
     */
    void dynamic_update(std::string &path);
};

#endif //TEMPUS_TEMPORAL_COMMUNITY_SEARCH_H
