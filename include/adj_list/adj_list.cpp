#include "adj_list.h"

#include <fstream>
#include <cinttypes>
#include <chrono>
#include <iostream>

AdjList::AdjList() : dynamicConnectivity(k){}

std::vector<std::set<int64_t>>
AdjList::computeComponents(std::unordered_set<int64_t> &uniqueVerticesSet) const{

    std::map<sequence::Id, std::set<int64_t>> componentsMap;
    for (const auto vertex : uniqueVerticesSet) {
        auto element = dynamicConnectivity.getId(vertex);
        componentsMap[element].emplace(vertex);
    }
    std::vector<std::set<int64_t>> components;
    components.reserve(componentsMap.size());
    for (const auto& [id, set] : componentsMap) {
        components.emplace_back(set);
    }
    return components;
}


std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>>
AdjList::addFromFile(const std::string &path) {
    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>> addMap, delMap;

    std::ifstream file(path);
    if (file.is_open()) {
        std::string command;
        int64_t source, destination, time;

        while (file >> command >> source >> destination >> time) {
            if (source == destination || source >= k || destination >= k) continue;
            if (command == "add") {
                addMap[time].emplace(source, destination);
            }
            /*
            if (command == "delete") {
                if (delMap.find(time) == delMap.end()) delMap.emplace(time, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>{});
                delMap.at(time).emplace(source, destination);
            }
             */
        }
        file.close();

        computeBatch(addMap, delMap);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time of addFromFile has taken " << duration.count() << "ms.\n" << std::endl;
    return addMap;
}

void AdjList::computeBatch(std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>> &addMap,
                           std::unordered_map<int64_t, std::unordered_set<UndirectedEdge, UndirectedEdgeHash>> &delMap) {

    for (const auto& [time, edgeData] : addMap) {
        for (const UndirectedEdge &edge : edgeData){
            if (!dynamicConnectivity.HasEdge(edge)){
                dynamicConnectivity.AddEdge(edge);
            }
            neighbourMap[edge.first][edge.second].emplace(time);
            neighbourMap[edge.second][edge.first].emplace(time);
            vertices.emplace(edge.first);
            vertices.emplace(edge.second);
        }
    }

    //TODO
    /*
    for (const auto& [time, edgeData] : delMap) {
        if (graph.find(time) == graph.end()) return;
        for (const UndirectedEdge &edge : edgeData){
            if (graph.at(time).HasEdge(edge)){
                graph.at(time).DeleteEdge(edge);
                edges.at(time).erase(edge);
            }
        }
    }
     */
}