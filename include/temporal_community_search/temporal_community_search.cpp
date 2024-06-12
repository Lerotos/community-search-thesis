#include "temporal_community_search.h"
#include <cinttypes>
#include <queue>
#include <algorithm>
#include <limits>
#include <iostream>

struct Timer{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::string fctName;
    explicit Timer(std::string& fctName){
        start = std::chrono::high_resolution_clock::now();
        this->fctName = fctName;
    }

    ~Timer(){
        end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Time of " << fctName << " has taken " << duration.count() << "ms.\n" << std::endl;
    }
};

tbb::global_control c(tbb::global_control::max_allowed_parallelism, 2);
temporal_community_search::temporal_community_search(std::string &path, int64_t theta, int64_t k, int64_t tau) {
    adjList.addFromFile(path);
    this->theta = theta;
    this->k = k;
    this->tau = tau;
}

void temporal_community_search::branchBound(std::set<int64_t> &S){
    auto uniqueVertices = adjList.vertices;

    if (uniqueVertices.size() < size) return;

    std::vector<std::set<int64_t>> components = adjList.computeComponents(uniqueVertices);

    for (auto& C : components) {
        //checks if S is a subset of C
        if (std::includes(C.begin(), C.end(), S.begin(), S.end())){

            if (C.size() < size) return;

            int64_t lc = meta_interval_intersection(C);

            if (lc >= tau){
                R = C;
                size = (int64_t) C.size();
                return;
            }

            int64_t ls = meta_interval_intersection(S);
            if (ls >= tau && S.size() < C.size()){
                int64_t u;

                //pick random node that is in C but not in S
                for (const int64_t node : C) {
                    if (S.count(node) == 0){
                        u = node;
                        break;
                    }
                }

                //first half
                int64_t fu = degreeCache[u];
                degreeCache[u] = 0;
                std::queue<int64_t> Q;
                Q.push(u);

                branchBoundHelper(Q, S);

                //second half
                degreeCache[u] = fu;

                S.emplace(u);
                int64_t lsTwo = meta_interval_intersection(S);
                if (lsTwo >= tau){
                    std::unordered_map<int64_t, int64_t> fCache;
                    std::queue<int64_t> B;
                    auto reducedSet = C;
                    for (const int64_t node : S) {
                        reducedSet.erase(node);
                    }
                    for (const int64_t v : reducedSet) {
                        auto tempS = S;
                        tempS.emplace(v);

                        int64_t lv = meta_interval_intersection(tempS);
                        if (lv < tau){
                            fCache[v] = degreeCache[v];
                            degreeCache[v] = 0;
                            B.push(v);
                        }
                    }
                    std::queue<int64_t> BQ = B;

                    branchBoundHelper(BQ, S);

                    while (!B.empty()){
                        degreeCache[B.front()] = fCache[B.front()];
                        B.pop();
                    }
                }
                S.erase(u);
            }
        }
    }
}

void temporal_community_search::branchBoundHelper(std::queue<int64_t> &Q, std::set<int64_t> &S){
    bool removeFlag = true;
    std::vector<int64_t> removeVector;
    while (!Q.empty()){
        int64_t u = Q.front();
        Q.pop();
        removeVector.emplace_back(u);
        auto data = removeNode(u, S);

        if (!data) {
            removeFlag = false;
            break;
        }

        auto H = data.value();
        for (const auto v : H) {
            Q.push(v);
        }
    }

    if (removeFlag) {
        //this deletes nodes from the graphs which is O(log^2 n) but I don't have a better way of doing this for now
        std::vector<UndirectedEdge> deletionCache = deleteNodes(removeVector);
        branchBound(S);
        //re-adding the removed nodes, similarly slow as removing them
        reAddNodes(deletionCache);
    }
}

template<typename T>
std::vector<int64_t>
temporal_community_search::genVertexVector(const T& Nu) const{
    std::vector<int64_t> vertexVector;
    for (const auto &[v, data] : Nu) {
        if (!deletedNodes.contains(v)){
            vertexVector.emplace_back(v);
        }
    }
    return vertexVector;
}

std::vector<UndirectedEdge> temporal_community_search::deleteNodes(std::vector<int64_t> &A){
    std::vector<UndirectedEdge> deletionCache;
    for (const auto u : A) {
        deleteSingleNode(u, deletionCache);
    }
    return deletionCache;
}

void temporal_community_search::deleteSingleNode(int64_t u, std::vector<UndirectedEdge> &deletionCache){
    auto vertexVector = genVertexVector(adjList.neighbourMap[u]);
    deletedNodes.emplace(u);
    adjList.vertices.erase(u);
    for (const int64_t v : vertexVector) {
        auto tempEdge = UndirectedEdge(u, v);
        if (adjList.dynamicConnectivity.HasEdge(tempEdge)){
            adjList.dynamicConnectivity.DeleteEdge(tempEdge);
            deletionCache.emplace_back(tempEdge);
        }
    }
}

void temporal_community_search::reAddNodes(std::vector<UndirectedEdge> &deletionCache){
    std::unordered_set<int64_t> tempDel;
    for (const auto& edge : deletionCache) {
        adjList.dynamicConnectivity.AddEdge(edge);
        if (deletedNodes.contains(edge.first)) {
            deletedNodes.erase(edge.first);
            tempDel.emplace(edge.first);
            adjList.vertices.emplace(edge.first);
        }
        if (deletedNodes.contains(edge.second)) {
            deletedNodes.erase(edge.second);
            tempDel.emplace(edge.second);
            adjList.vertices.emplace(edge.second);
        }
    }
    for (auto n:tempDel) {
        addNode(n);
    }
}


std::optional<std::unordered_set<int64_t>>
temporal_community_search::removeNode(int64_t u, std::set<int64_t> &S){
    bool flag = true;
    std::unordered_set<int64_t> resultSet;
    tbb::spin_mutex mtx;

    auto f = [&mtx, &resultSet, &S, &flag](int64_t v){
        tbb::spin_mutex::scoped_lock l(mtx);
        resultSet.emplace(v);
        if (S.count(v)){
            flag = false;
        }
    };

    removeNodeHelper(u, f);
    if (!flag) return {};
    return resultSet;
}

template <typename F>
void temporal_community_search::removeNodeHelper(int64_t u, F&& f){
    auto Nu = adjList.neighbourMap[u];
    std::vector<int64_t> vertexVector = genVertexVector(Nu);
    auto func = [&](const tbb::blocked_range<size_t>& r){
        for (auto i = r.begin(); i < r.end(); ++i) {
            int64_t v = vertexVector[i];
            auto times = Nu[v];

            auto MIv = MICache[v];

            for (const auto t : times) {
                std::map<int64_t, int64_t>::iterator itStart, itEnd;
                itStart = MIv.lower_bound(t);
                itEnd = MIv.lower_bound(t + theta);

                int64_t j = std::distance(MIv.begin(), itStart);

                for (auto it = itStart; it != itEnd; ++it){
                    DCache[v][j]--;
                    auto DC = DCache[v][j];
                    if (DC + 1 >= k && DC < k){
                        int64_t temp = degreeCache[v];
                        int64_t ts = it->first;
                        int64_t te = it->second;
                        degreeCache[v] -= (te - ts);
                        if (degreeCache[v] == theta) degreeCache[v] = 0;

                        if (temp >= tau && degreeCache[v] < tau){
                            f(v);
                        }
                    }
                    j++;
                }
            }
        }
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexVector.size()), func);
}

void temporal_community_search::addNode(int64_t u){
    auto Nu = adjList.neighbourMap[u];
    std::vector<int64_t> vertexVector = genVertexVector(Nu);

    auto func = [&](const tbb::blocked_range<size_t>& r){
        for (auto i = r.begin(); i < r.end(); ++i) {

            int64_t v = vertexVector[i];
            if (preCutNodes.contains(v)) continue;
            auto times = Nu[v];

            auto MIv = MICache[v];
            auto Dv = DCache[v];

            for (const auto t : times) {
                std::map<int64_t, int64_t>::iterator itStart, itEnd;
                itStart = MIv.lower_bound(t);
                itEnd = MIv.lower_bound(t + theta);

                int64_t j = std::distance(MIv.begin(), itStart);

                for (auto it = itStart; it != itEnd; ++it) {
                    DCache[v][j]++;
                    if (DCache[v][j] - 1 < k && DCache[v][j] >= k){
                        int64_t ts = it->first;
                        int64_t te = it->second;
                        if (degreeCache[v] == 0) degreeCache[v] = theta;
                        degreeCache[v] += (te - ts);
                    }
                    j++;
                }
            }
        }
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexVector.size()), func);
}

int64_t temporal_community_search::meta_interval_intersection(std::set<int64_t> &components){

    std::map<int64_t, int64_t> MIR;
    MIR.emplace(std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max());
    int64_t len = std::numeric_limits<int64_t>::max() - theta;

    for (const int64_t u : components) {
        std::map<int64_t, int64_t> MIT;
        int64_t tlen = 0;
        auto MIu = MICache[u];
        auto Du = DCache[u];
        auto p = (int64_t) MIR.size();
        int64_t j = 1;
        int64_t i = 0;
        for (const auto& [startMIu, endMIu] : MIu) {
            if (Du[i] >= k){


                for ( const auto& [startMIR, endMIR] : MIR) {
                    if (j > p) break;

                    //first while loop
                    if (endMIR <= startMIu){
                        j++;
                        continue;
                    }

                    //second while loop
                    if (endMIR <= endMIu) {
                        int64_t l = std::max(startMIR, startMIu);
                        int64_t r = endMIR;
                        if (l < r){
                            MIT.emplace(l, r);
                            tlen += (r-l);
                        }
                        j++;
                        continue;
                    }

                    //if condition after both while loops
                    if (startMIR <= endMIu){
                        int64_t l = std::max(startMIR, startMIu);
                        int64_t r = endMIu;
                        if (l < r){
                            MIT.emplace(l, r);
                            tlen += (r-l);
                        }
                    }
                    break;

                }

            }
            i++;
        }
        MIR.clear();
        MIR.insert(MIT.begin(), MIT.end());
        len = tlen;
    }

    if (len > 0) len += theta;

    return len;

}

void temporal_community_search::meta_interval_decomposition(int64_t u){

    std::map<int64_t, int64_t> T;
    auto Nu = adjList.neighbourMap[u];
    for (const auto& [neighbour, times] : Nu) {

        int64_t carryoverTime = std::numeric_limits<int64_t>::min();
        for (const int64_t time: times) {
            //handles repeated neighbours
            if (time < carryoverTime){
                carryoverTime = time;
            }

            //inserts {t + theta, -1} in the following iteration because of repeated neighbours
            if (carryoverTime > std::numeric_limits<int64_t>::min()){
                if(T.find(carryoverTime) == T.end()) T.emplace(carryoverTime, 0);
                T[carryoverTime]--;
            }

            carryoverTime = time + theta;

            if(T.find(time) == T.end()) T.emplace(time, 0);
            T[time]++;
        }
        //inserts {t + theta, -1} for the last iteration
        if(T.find(carryoverTime) == T.end()) T.emplace(carryoverTime, 0);
        T[carryoverTime]--;
    }

    //          sum
    std::vector<int64_t> D;

    //       start    end
    std::map<int64_t, int64_t> MI;
    int64_t carryOver = -1;
    int64_t d = 0;
    for (const auto [time, sum] : T) {
        //populate MI
        if (d > 0 && carryOver > -1) MI.emplace(carryOver, time);

        //populate D
        d += sum;
        if (d > 0) D.emplace_back(d);

        carryOver = time;

    }

    MICache[u] = MI;
    DCache[u] = D;
}


void temporal_community_search::TGR(){
    std::string s = "TGR";
    Timer timer(s);

    preCutNodes.clear();

    std::queue<int64_t> Q{};
    if (adjList.neighbourMap.empty()) return;
    auto vertexVector = genVertexVector(adjList.neighbourMap);
    tbb::spin_mutex mtx;

    auto func = [&](const tbb::blocked_range<size_t>& r){
        for (auto i = r.begin(); i < r.end(); ++i) {
            int64_t u = vertexVector[i];
            auto Nu = adjList.neighbourMap[u];

            meta_interval_decomposition(u);
            int64_t degreePersistence = computeDegreePersistence(u);
            degreeCache[u] = degreePersistence;

            if (degreePersistence < tau) {
                tbb::spin_mutex::scoped_lock l(mtx);
                Q.push(u);

                deleteSingleNode(u, TGRDeletionCache);
                preCutNodes.emplace(u);
            }
        }
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexVector.size()), func);

    tbb::concurrent_unordered_set<int64_t> tn;

    auto f = [&](int64_t v) {
        Q.push(v);
        tn.emplace(v);
    };

    while (!Q.empty()){
        int64_t u = Q.front();
        Q.pop();
        removeNodeHelper(u, f);
    }

    for (auto v:tn) {
        deleteSingleNode(v, TGRDeletionCache);
    }
}

int64_t temporal_community_search::computeDegreePersistence(int64_t u){
    auto MI = MICache[u];
    auto D = DCache[u];
    int64_t degreePersistence = 0;
    int64_t i = 0;
    for (const auto& [start, end] : MI) {
        if (D[i] >= k) degreePersistence += (end - start);
        i++;
    }
    if (degreePersistence > 0) {
        degreePersistence += theta;
    }
    return degreePersistence;
}

void temporal_community_search::reset(){
    reAddNodes(TGRDeletionCache);
    TGRDeletionCache.clear();

    size = 0;
    R.clear();
}

void temporal_community_search::meta_interval_decomposition_update(int64_t u, int64_t time, int64_t neighbour,
                                                                   std::set<int64_t> &neighbourTimes) {

    //used to later add theta to the degree if it hasn't happened yet
    bool zeroFlag = false;
    if (degreeCache[u] == 0) zeroFlag = true;

    bool startFlag, endFlag;
    startFlag = endFlag = true;
    int64_t timeEnd = time + theta;
    auto MIBegin = MICache[u].begin();
    auto MIEnd = MICache[u].end();

    //Checking if new entries need to be generated at the start and end time
    if (MICache[u].count(time))
        startFlag = false;
    if (MICache[u].count(timeEnd))
        endFlag = false;

    int64_t dEnd = 0, dStart = 0, d = 0;

    //setting up the iterators
    std::map<int64_t, int64_t>::iterator itStart, itEnd, itCache;
    itStart = MICache[u].lower_bound(time);
    itEnd = MICache[u].upper_bound(timeEnd);
    itCache = MIEnd;

    //catches edges that start and end at the same time as a previous edge
    if (MICache[u].lower_bound(time)->second == timeEnd && !startFlag){
        endFlag = false;
        dEnd++;
        degreeCache[u] += theta;
    }

    //sets up itStart before the insert spot if possible, adjusts dStart to not be affected by it
    if (itStart != MIBegin && itStart != MIEnd ){//&& itStart->first != time
        --itStart;
        dStart++;
    }

    //setting up boundaries for DCache updates
    dStart += std::distance(MIBegin, itStart);
    dEnd += std::distance(MIBegin, itEnd) + 2;

    bool verySpecificCornerCaseFlag = false;
    bool noRepeatingNeighbourAfter = true;
    if (!neighbourTimes.empty()){
        auto it = neighbourTimes.lower_bound(time - theta);
        int64_t end, start = -1;
        while (it != neighbourTimes.end() && *it < time + theta){
            if (start == -1) start = *it;
            end = *it + theta;
            ++it;
        }
        if (start != -1){
            //identical node, can be discarded
            if (start == time) return;

            //existing node comes after new entry, reduce end time of new entry
            if (end > time + theta) {
                timeEnd = MICache[u].lower_bound(time)->first;
                endFlag = false;
                if (neighbourTimes.count(timeEnd) && MICache[u].upper_bound(timeEnd) != itEnd) dEnd--;
                verySpecificCornerCaseFlag = true;
            }

            //existing node comes before new entry
            //reduce end time of existing node
            //adjust d
            if (start < time){
                auto tempIt = --MICache[u].upper_bound(time);
                auto timesIt = --neighbourTimes.lower_bound(time);

                //case: non-repeating neighbour starts or ends during the interval
                if (DCache[u][dStart - 1] > 1){
                    auto temp = tempIt->second;
                    tempIt->second = time;
                    MICache[u][time] = temp;

                    startFlag = false;
                    d = DCache[u][dStart - 1];
                    DCache[u].insert(DCache[u].begin() + dStart, d);

                    while (tempIt->second != *timesIt + theta){
                        ++tempIt;
                    }
                    auto cache = tempIt;
                    ++tempIt;

                    //case: last interval before entry is from a repeated neighbour
                    if (tempIt->first != cache->second){
                        --tempIt;
                        --tempIt;
                        cache->second = tempIt->second;
                        noRepeatingNeighbourAfter = false;
                        itEnd = --MICache[u].upper_bound(timeEnd);
                        itStart = MICache[u].lower_bound(temp);

                    //case: last interval before entry is from a non-repeating neighbour
                    }else{
                        dStart = std::distance(MIBegin, tempIt);

                        cache->second = tempIt->second;
                        MICache[u].erase(tempIt->first);
                        DCache[u].erase(DCache[u].begin() + dStart);
                        itEnd = MICache[u].upper_bound(timeEnd);
                        itStart = MICache[u].lower_bound(time);
                    }

                    dEnd = std::distance(MIBegin, itEnd) + 2;
                    dStart = std::distance(MIBegin, MICache[u].lower_bound(temp));

                    //reduces the degree by the sub-interval cut from the repeating neighbour
                    if (d == k){
                        int64_t dif = temp - cache->first;
                        degreeCache[u] -= dif;
                    }
                } else {    //case: no non-repeating neighbours start ot end during the interval
                    tempIt->second = time;
                }
            }
        }
    }
    neighbourTimes.emplace(time);

    //reducing dEnd in corner cases
    if (!endFlag){
        dEnd--;
        auto temp = itEnd;
        --temp;
        if (temp->second != itEnd->first && startFlag) verySpecificCornerCaseFlag = true;
        if (temp->second == itEnd->first || verySpecificCornerCaseFlag){
            dEnd--;
        }
    }
    if (!startFlag){
        dEnd--;
        if (!endFlag) dEnd--;
    }

    int64_t controlD = -1;
    //updating dCache
    //updates d to the value of the previous entry if it was recent enough
    if (itStart->second > time && itStart->first <= time && itStart != itEnd){
        d = DCache[u][dStart - 1];
        controlD = d;
    }
    for (int64_t j = dStart; j < dEnd; ++j) {
        if (j == dStart && startFlag){  //inserts the first entry
            DCache[u].insert(DCache[u].begin() + dStart, ++d);
        } else if (j + 1 == dEnd && endFlag){ //inserts the last entry
            if (d > 1 && noRepeatingNeighbourAfter) {
                //catches the case of two edges ending at the same time within the interval
                if (controlD > 1){
                    if (itStart != MIBegin){
                        auto temp = itStart;
                        --temp;
                        if ((itStart->first == temp->second)){
                            if (DCache[u][dStart - 1] > DCache[u][dStart - 2] + 1)
                                d -= d - DCache[u][dStart - 2];
                        }
                        else d = 2;
                    }
                }

                DCache[u].insert(DCache[u].begin() + j, d - 1);
            }
        } else {    //increases all values inbetween
            DCache[u][j]++;
            d = DCache[u][j];
        }
    }

    int64_t i = std::max((int64_t) 0, dStart - 1);

    //used to fill gaps that are created through an insert
    std::unordered_map <int64_t, int64_t> gapMap;
    std::set<int64_t> gapSet;

    std::map<int64_t, int64_t>::iterator it;

    //inserts the start of the interval
    auto insert_start_func = [&](){
        if (it->first <= timeEnd && it != MIEnd){ //case: current it doesn't start after timeEnd
            MICache[u][time] = it->first;
        } else if (itCache->second <= time || itCache == itEnd){ //case: current it starts after timeEnd
            MICache[u][time] = timeEnd;
            endFlag = false;
        }

        int64_t subtraction = 0;
        //case: insert happens after the previous entry starts but before it ends
        if (itCache->second > time && itCache != MIEnd) {
            if (DCache[u][i - 1] >= k && DCache[u][i] >= k) {    //temporarily reduces the degree, this will later be added back
                subtraction = itCache->second - time;
            }

            int64_t temp = itCache->second;
            itCache->second = time;
            MICache[u][time] = temp;
        }
        startFlag = false;

        if (DCache[u][i] >= k){
            degreeCache[u] += (MICache[u][time] - time) - subtraction;
        }
    };

    //inserts the end of the interval
    auto insert_end_func = [&](){
        auto tempIt = it;
        --tempIt;

        if (tempIt->second < timeEnd){  //case: interval ends outside an existing interval
            MICache[u][tempIt->second] = timeEnd;
        } else {    //case: interval end inside an existing interval
            int64_t temp = tempIt->second;
            MICache[u][tempIt->first] = timeEnd;
            MICache[u][timeEnd] = temp;
        }
    };

    //increases the degree after an interval was changed
    auto degree_func = [&](){
        int64_t DC = DCache[u][i - 1];
        if (DC >= k && DC - 1 < k && itCache->second <= time + theta && itCache->first >= time){
            degreeCache[u] += (itCache->second - itCache->first);
        }
    };

    //updating MICache
    for (it = itStart; it != itEnd; ++it){
        if (startFlag && it->first > time){ //inserts before current it
            insert_start_func();
            i++; //moves it forward one to skip the new entry
        }

        if (endFlag && it->first > timeEnd ){
            insert_end_func();

            endFlag = false;
            break;
        }

        if (it != itStart){
            degree_func();

            auto tempIt = it;
            --tempIt;

            //used to fill intervals that were outside existing intervals before and are not adjacent to time or timeEnd
            if (tempIt->second != it->first && it->first <= timeEnd && tempIt->second >= time){
                gapMap.emplace(tempIt->second, it->first);
                gapSet.emplace(i);
            }
        }

        itCache = it;
        i++;
    }

    //last iteration to emulate it <= itEnd
    if (startFlag){
        insert_start_func();
    }
    if (endFlag){
        insert_end_func();
        degree_func();
    }

    //filling the gaps here to not mess with the iterators
    if (!gapMap.empty()){
        for (const auto [start, end]:gapMap) {
            MICache[u][start] = end;
        }
        for (const auto x : gapSet){
            DCache[u].insert(DCache[u].begin() + x, 1);
        }
    }

    //adds the baseline degree if this node didn't have a degree before
    if (zeroFlag && degreeCache[u] > 0) degreeCache[u] += theta;
}

void temporal_community_search::dynamic_update(std::string& path){
    std::string s = "dynamic_update";
    Timer timer(s);
    std::queue<int64_t> Q{};

    auto tempDel = preCutNodes;

    preCutNodes.clear();

    size = 0;
    R.clear();

    tbb::concurrent_unordered_set<int64_t> nodes;
    tbb::concurrent_unordered_set<int64_t> fromScratch;

    auto neighbourCache = adjList.neighbourMap;

    auto inserts = adjList.addFromFile(path);

    //restructures the data so it can be used in parallel
    std::unordered_map<int64_t, std::unordered_map<int64_t, std::set<int64_t>>> map;
    for (const auto& [time, edgeSet]:inserts) {
        for (const auto& edge: edgeSet) {
            map[edge.first][edge.second].emplace(time);
            map[edge.second][edge.first].emplace(time);
        }
    }

    auto vertexVector = genVertexVector(map);
    auto func = [&](const tbb::blocked_range<size_t>& r){
        for (auto i = r.begin(); i < r.end(); ++i) {

            int64_t v = vertexVector[i];
            auto innerMap = map[v];

            auto neighbourList = neighbourCache[v];

            for (const auto& [neighbour, times]:innerMap) {
                if (degreeCache.contains(v)){
                    for (const auto time:times) {
                        meta_interval_decomposition_update(v, time, neighbour, neighbourList[neighbour]);
                    }
                } else {
                    fromScratch.emplace(v);
                }
                nodes.emplace(v);
            }
        }
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexVector.size()), func);

    for (auto u : fromScratch){
        meta_interval_decomposition(u);
        int64_t degreePersistence = computeDegreePersistence(u);
        degreeCache.emplace(u, degreePersistence);
    }

    for (auto node:tempDel) {
        nodes.emplace(node);
    }

    for (const auto n : nodes) {
        if (degreeCache[n] < tau) {
            deleteSingleNode(n, TGRDeletionCache);
            preCutNodes.emplace(n);
            Q.push(n);
        }
    }

    tbb::spin_mutex m;
    tbb::concurrent_unordered_set<int64_t> tn;
    auto f = [&](int64_t v) {
        Q.push(v);
        tn.emplace(v);
    };

    while (!Q.empty()){
        int64_t u = Q.front();
        Q.pop();
        removeNodeHelper(u, f);
    }

    for (auto v:tn) {
        deleteSingleNode(v, TGRDeletionCache);
    }
}
