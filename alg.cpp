#include "stdafx.h"

// code from OPIM
double Alg::MaxCoverVanilla(const int targetSize)
{
    // optimization with minimum upper bound among all rounds [Default].
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    // degMap: map degree to the nodes with this degree
    RRsets degMap(maxDeg + 1);

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// In the case of identical marginal，the node with largest out degree is chosen
double Alg::MaxCoverOutDegPrority(const int targetSize)
{
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// max cover used in the IM-sentinel Phase
double Alg::MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize)
{
    // seedSet: the sentinel set obtained in the Sentinel Set Selection Phase
    std::unordered_set<uint32_t> subSeedSet(seedSet.begin(), seedSet.end());
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto node : seedSet)
    {
        sumInf += coverage[node];
        _vecSeed.push_back(node);
        coverage[node] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[node])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                coverage[nodeIdx]--;
            }
        }
    }

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverSentinelSet(const int targetSize, const int totalTargetSize)
{
    //targetSize: the size of the sentinel set
    //totalTargetSize: the total number of the seed set
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = totalTargetSize;
                auto degBound = deg;
                FRset vecBound(totalTargetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                goto afterGreedy;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            _vecVldtInf.push_back(sumInf);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

afterGreedy:
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    // std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", seed set: " << _vecSeed.size() << ", min-bound: " << _boundMin <<
    //           ", last-bound: " << _boundLast << std::endl;

    if (_vecSeed.size() == targetSize)
    {
        // if the sample size is sufficiently large, the result is reliable
        if (_numRRsets > 1000)
        {
            return finalInf;
        }
    }

    // if covering all the RR-sets, the sentinel set may include some nodes which cover only a small number of RR-sets.
    // such nodes should not be included.
    uint32_t threshold = 0.9 * sumInf;
    for (int i = _vecSeed.size() - 1; i > 0; i--)
    {
        if (_vecVldtInf[i - 1] >= threshold)
        {
            degMap[0].push_back(_vecSeed[i]);
            _vecSeed.pop_back();
        }
        else
        {
            break;
        }
    }

    // std::cout << "seedset size reaching 0.9 coverage: " << _vecSeed.size() << std::endl;
    // the following code is to select the nodes with large out-degree. 
    int seedSetSize = (degMap[0].size() > targetSize) ? targetSize : degMap[0].size();
    std::vector<std::pair<uint32_t, uint32_t>> vecHeap;

    for (int i = 0; i < seedSetSize; i++)
    {
        auto &node = degMap[0][i];
        vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
    }

    std::make_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    const auto nodeNum = degMap[0].size();


    for (int i = seedSetSize; i < nodeNum; i++)
    {
        uint32_t node = degMap[0][i];
        uint32_t currDeg = _vecOutDegree[node];

        if (currDeg > vecHeap[0].first)
        {
            std::pop_heap(vecHeap.begin(), vecHeap.end());
            vecHeap.pop_back();
            vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
            std::push_heap(vecHeap.begin(), vecHeap.end());
        }
    }

    std::sort_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    std::unordered_set<uint32_t> seedHashSet(_vecSeed.begin(), _vecSeed.end());

    for (auto &node : vecHeap)
    {
        if (_vecSeed.size() >= targetSize)
        {
            break;
        }

        if (seedHashSet.find(node.second) != seedHashSet.end())
        {
            std::cout << "node exist" << std::endl;
            continue;
        }

        _vecSeed.push_back(node.second);
        seedHashSet.insert(node.second);

        if (_vecVldtInf.size() < _vecSeed.size())
        {
            _vecVldtInf.push_back(sumInf);
        }
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverTopK(const int targetSize)
{
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        //if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }

    Nodelist sortedNode(_numV); // sortedNode: record the sorted nodes in ascending order of degree
    Nodelist nodePosition(_numV); // nodePosition: record the position of each node in the sortedNode
    Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;

    for (auto& nodes : degMap)
    {
        degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;

        for (auto& node : nodes)
        {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }

    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;

    for (auto deg = maxDeg + 1; deg--;)
    {
        if (degreePosition[deg] <= _numV - targetSize)
        {
            sumTopk += deg * (degreePosition[deg + 1] - (_numV - targetSize));
            break;
        }

        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }

    _boundMin = 1.0 * sumTopk;
    _vecSeed.clear();
    size_t sumInf = 0;

    /*
    * sortedNode: position -> node
    * nodePosition: node -> position
    * degreePosition: degree -> position (start position of this degree)
    * coverage: node -> degree
    * e.g., swap the position of a node with the start position of its degree
    * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
    */
    for (auto k = targetSize; k--;)
    {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        _vecSeed.push_back(seed);
        coverage[seed] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[seed])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                const auto currPos = nodePosition[nodeIdx]; // The current position
                const auto currDeg = coverage[nodeIdx]; // The current degree
                const auto startPos = degreePosition[currDeg]; // The start position of this degree
                const auto startNode = sortedNode[startPos]; // The node with the start position
                // Swap this node to the start position with the same degree, and update their positions in nodePosition
                std::swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;

                // If the start position of this degree is in top-k, reduce topk by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }

        _boundLast = 1.0 * (sumInf + sumTopk);

        if (_boundMin > _boundLast) _boundMin = _boundLast;
    }

    _boundMin *= 1.0 * _numV / _numRRsets;
    _boundLast *= 1.0 * _numV / _numRRsets;
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << _boundMin <<
              ", last-bound: " << _boundLast << '\n';
    return finalInf;
}

double Alg::MaxCover(const int targetSize)
{
    if (targetSize >= 1000) return MaxCoverTopK(targetSize);

    return MaxCoverVanilla(targetSize);
}

double Alg::MaxWCover(const int targetSize) {
    // for weighted max cover
    std::priority_queue<std::pair<uint32_t, double>, std::vector<std::pair<uint32_t, double>>, CompareBySecondDouble> heap;
    std::vector<double> coverage(_numV, 0.0);
    for (uint32_t i=0; i<_numV; i++) {
        coverage[i] = this->_hyperGraph._nodeGain[i];
        heap.push(std::make_pair(i, coverage[i]));
    }

    _vecSeed.clear();
    uint32_t max_idx = 0;
    double cov_num = 0.0;
    std::vector<bool> vecBoolVst = std::vector<bool>(this->_numRRsets, false);
    while (_vecSeed.size() < targetSize)
    {
        std::pair<uint32_t, double> top = heap.top();
        heap.pop();
        
        // Lazy Update
        if (top.second > coverage[top.first]) {
            // Update coverage of top
            top.second = coverage[top.first];
            heap.push(top);
            continue;
        }
        
        max_idx = top.first;
        cov_num += coverage[max_idx];
        _vecSeed.push_back(max_idx);
        // LogInfo("node", max_idx);
        // LogInfo("coverage", coverage[max_idx]);
        
        double origin_cov = coverage[max_idx];

        std::vector<size_t>& RRsets_cov_by_max_node = this->_hyperGraph._FRsets[max_idx];
        // LogInfo("Updating");
        for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
            size_t RRset_idx = RRsets_cov_by_max_node[j];
            if (vecBoolVst[RRset_idx]) {
                continue;
            }
            vecBoolVst[RRset_idx] = true;
            Nodelist& node_list = this->_hyperGraph._RRsets[RRset_idx];
            // timer.get_operation_time();
            for (uint32_t& node: node_list){
                if (node != max_idx) {
                    coverage[node] -= this->_hyperGraph._RRweights[RRset_idx];
                }
            }
        }

        //update
        coverage[max_idx] = 0.0;
    }
    // σ̂(S) = (n + W_out) × Σ_{covered r} w_r / R  (importance-corrected; _norm_const per wrr mode)
    return cov_num * _norm_const / _numRRsets;
}
void Alg::set_prob_dist(const ProbDist dist)
{
    _probDist = dist;
    _hyperGraph.set_prob_dist(dist);
    _hyperGraphVldt.set_prob_dist(dist);
}

void Alg::set_vanilla_sample(const bool isVanilla)
{
    if (isVanilla)
    {
        std::cout << "Vanilla sampling method is used" << std::endl;
    }

    _hyperGraph.set_vanilla_sample(isVanilla);
    _hyperGraphVldt.set_vanilla_sample(isVanilla);
}

void Alg::set_wrr_sample_mode(const WRRSampleMode mode)
{
    _hyperGraph.set_wrr_sample_mode(mode);
    _hyperGraphVldt.set_wrr_sample_mode(mode);
    _norm_const = static_cast<size_t>(_hyperGraph.get_wrr_norm_const());
    switch (mode) {
    case WRR_WINDEG:
        std::cout << "WRR root sampling: weighted in-degree" << std::endl;
        break;
    case WRR_WOUTDEG:
        std::cout << "WRR root sampling: weighted out-degree (+1)" << std::endl;
        break;
    case WRR_OUTDEG:
        std::cout << "WRR root sampling: out-degree (+1)" << std::endl;
        break;
    case WRR_UNIFORM:
        std::cout << "WRR root sampling: uniform" << std::endl;
        break;
    case WRR_INDEG:
    default:
        std::cout << "WRR root sampling: in-degree" << std::endl;
        break;
    }
}

double Alg::EfficInfVldtAlg()
{
    return EfficInfVldtAlg(_vecSeed);
}

double Alg::EfficInfVldtAlg(const Nodelist vecSeed)
{
    Timer EvalTimer("Inf. Eval.");
    std::cout << "  >>>Evaluating influence in [0.99,1.01]*EPT with prob. 99.9% ...\n";
    const auto inf = _hyperGraph.EfficInfVldtAlg(vecSeed);
    //const auto inf = _hyperGraphVldt.EfficInfVldtAlg(vecSeed);
    std::cout << "  >>>Down! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
    return inf;
}

double Alg::estimateRRSize()
{
    const int sampleNum = 100;
    _hyperGraph.BuildRRsets(sampleNum);
    double avg= _hyperGraph.HyperedgeAvg();
    _hyperGraph.RefreshHypergraph();
    return avg;
}

double Alg::subsimOnly(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _hyperGraphVldt.BuildRRsets(numR); // R2
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCover(targetSize);
        time2 += timerSubsim.get_operation_time();
        const auto infVldt = _hyperGraphVldt.CalculateInf(_vecSeed);

        const auto degVldt = infVldt * _numRRsets / _numV;
        auto upperBound = _boundMin;
        

        // const auto upperDegOPT = upperBound * _numRRsets / _numV;
        auto upperDegOPT = (double)degVldt / (1-1/e);
        const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double avgSize = _hyperGraph.HyperedgeAvg();

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            return 0;
        }
    }

    return 0.0;
}

double Alg::subsimWeight(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM-W");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildWRRsets(numR); // R1
        _hyperGraphVldt.BuildWRRsets(numR); // R2
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxWCover(targetSize);
        time2 += timerSubsim.get_operation_time();
        const auto infVldt = _hyperGraphVldt.CalculateInf_W(_vecSeed);
        LogInfo("infSelf", infSelf);
        LogInfo("infVldt", infVldt);
        const auto degVldt = infVldt * _numRRsets / _norm_const;
        auto upperBound = _boundMin;
        auto degSelf = infSelf * _numRRsets / _norm_const;
        
        const auto upperDegOPT = upperBound * _numRRsets /_norm_const;
        const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        const auto upperOPT = pow2(sqrt(degSelf / (double)(1-1/e) + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _norm_const / _numRRsets) << ", upperBound: " << (upperOPT * _norm_const / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double avgSize = _hyperGraph.HyperedgeAvg();

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            return 0;
        }
    }

    return 0.0;
}

int decideMultiple(int ratio, int numRRsets)
{
    int multiple = 1;

    if (numRRsets < 100)
    {
        return multiple;
    }

    if (ratio >= 32)
    {
        multiple = 8;
    }
    else if (ratio >= 16)
    {
        multiple = 4;
    }
    else if (ratio >= 4)
    {
        multiple = 2;
    }
    else
    {
        multiple = 1;
    }

    return multiple;
}

double Alg::subsimWithTrunc(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(3 * log(1 / delta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    int multiple = 1;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCoverOutDegPrority(targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();
        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR sets: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR set: " << time4 << std::endl;
            std::cout << "==>Time for greedy: " << time2 << std::endl;
            return 0;
        }
    }

    return 0.0;
}

double Alg::IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr)
{
    size_t vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    size_t R1RRsets = _hyperGraph.get_RR_sets_size();
    int multiple = 4;
    double estimateAppr = 0.0;
    double lowerSelect = 0;
    int maxMultiple = 3;
    double infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    double degVldt = infVldt * vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt * multiple + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0) ;
    estimateAppr = (lowerSelect / (multiple * vldtRRsets)) / (upperOPT / R1RRsets);

    if (estimateAppr < targetAppr)
    {
        return 0.0;
    }

    _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, vldtRRsets * multiple);
    infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    degVldt = infVldt *  vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0);
    double newAppr = (lowerSelect / vldtRRsets) / (upperOPT / R1RRsets);
    return (newAppr > targetAppr) ? newAppr : 0.0;
}

double Alg::FindRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    size_t subSeedSetSize = _vecSeed.size();
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    // delta for upper bound on the number of RR sets
    const double delta_upper = delta / 3.0;
    const double alpha = sqrt(log(3.0 / delta_upper));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize - subSeedSetSize) + log(3.0 / delta_upper)));
    //const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto numRbase = size_t(_baseNumRRsets);
    const auto maxNumR = size_t(2.0 * _numV * pow2(alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    double currApprox = 0.0;
    double infSelf = 0.0;
    int multiple = 1;
    std::unordered_set<uint32_t> subSeedSet(_vecSeed.begin(), _vecSeed.end());
    std::vector<uint32_t> vecSubSeed(_vecSeed.begin(), _vecSeed.end());

    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsetsEarlyStop(subSeedSet, numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        infSelf = MaxCoverIMSentinel(vecSubSeed, targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();

        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - targetEpsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR in IM-Sentinel phase: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR in IM-Sentinel phase: " << time4 << std::endl;
            std::cout << "==>Time for greedy in IM-Sentinel phase: " << time2 << std::endl;
            std::cout << "==>Influence via R2 in IM-Sentinel phase: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            return 0;
        }
    }

    return 0.0;
}

double Alg::FindDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double x = (1.0 - 1.0 / totalTargetSize);
    const int minSubSize = ceil(log(1 - epsilon) / log(x));
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
    const auto numRbasePrevious = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize);
    const auto numRbase = size_t(_baseNumRRsets);
    
    // the successful probability of at least 1-delta/3
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    int multiple = 1;
    double infVldt = 0.0;
    bool firstRound = true;


    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        //build R1
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        _vecVldtInf.clear();
        int targetSize = firstRound ? (totalTargetSize / 4) : (totalTargetSize / 8);
        firstRound = false;
        const auto infSelf = MaxCoverSentinelSet(targetSize, totalTargetSize);
        time2 += timerSubsim.get_operation_time();
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize)
        {
            //low influence
            continue;
        }

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            infVldt = _vecVldtInf[i];
            double lowerDeg = infVldt;
            double upperDeg = _boundMin;
            double a = log(numIter * 6.0  / delta);
            double lower = pow2(sqrt(lowerDeg + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
            double upper = pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            upper = (upper > numR) ? numR : upper;
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2;

            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }


        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        _hyperGraphVldt.RefreshHypergraph();
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        time4 += timerSubsim.get_operation_time();
        double degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        double upperDegOPT = upperBound * _numRRsets / _numV;
        double lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;

        if (lowerSelect < 0)
        {
            lowerSelect = 1.0 * _vecSeed.size() / _numV * _numRRsets * multiple;
        }

        double upperOPT = pow2(sqrt(upperDegOPT + a1 / 2.0) + sqrt(a1 / 2.0));
        upperOPT = (upperOPT > _numRRsets) ? _numRRsets : upperOPT;
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / (_numRRsets)) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr)
        {
            goto succ;
        }

        if (_numRRsets < 100)
        {
            continue;
        }

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();

        if (fullRRSize / truncRRSize < 2)
        {
            continue;
        }

        double lowerThreshold = (upperOPT * _numV / _numRRsets) * targetAppr;

        if ((1.0 * infVldt / multiple) > lowerThreshold && lowerThreshold > 0)
        {
            double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);
            time4 += timerSubsim.get_operation_time();

            if (newAppr > targetAppr)
            {
                std::cout << "increase R2 successfully" << std::endl;
                infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
                goto succ;
            }
        }
    }

succ:
    std::cout << "==>Time for full RR in SentinelSet phase: " << time1  << std::endl;
    std::cout << "==>Time for truncated RR in SentinelSet phase: " << time4 << std::endl;
    std::cout << "==>Time for greedy in SentinelSet phase: " << time2 << std::endl;
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", inf: " << infVldt << std::endl;
    std::cout << "==>total time for SentinelSet phase: " << timerSubsim.get_total_time() << std::endl;
    return 0.0;
}

double Alg::subsimWithHIST(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    _baseNumRRsets = 3 * log(1 / delta);
    
    std::cout << std::endl;
    std::cout << "Sentinel Set Selection Phase" << std::endl;
    FindDynamSub(targetSize, epsilon / 2, delta / 2);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();

    std::cout << std::endl;
    std::cout << "IM-Sentinel Phase" << std::endl;
    FindRemSet(targetSize, epsilon / 2, epsilon, delta / 2);
    _res.set_running_time(timerSubsim.get_total_time());
    return 0.0;
}

// Fixed number of RR sets sampling and then run greedy algorithm
double Alg::fixed_subsim(const int targetSize, int num_samples)
{
    double time1 = 0.0, time2 = 0.0;
    Timer timer("fixed_subsim");
    
    // Build fixed number of RR sets
    _hyperGraph.BuildRRsets(num_samples);
    _hyperGraphVldt.BuildRRsets(num_samples);
    _numRRsets = _hyperGraph.get_RR_sets_size();
    time1 += timer.get_operation_time();
    
    // Run greedy algorithm (MaxCover)
    const auto infSelf = MaxCover(targetSize);
    time2 += timer.get_operation_time();

    // Set results (seed selection only)
    _res.set_running_time(timer.get_total_time());
    _res.set_influence_original(infSelf);
    _res.set_seed_vec(_vecSeed);
    _res.set_RR_sets_size(_numRRsets);
    _res.set_sampling_time(time1);
    _res.set_selection_time(time2);

    std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
    return 0.0;
}

// Fixed number of weighted RR sets sampling and then run weighted greedy algorithm
double Alg::fixed_subsimW(const int targetSize, int num_samples)
{
    double time1 = 0.0, time2 = 0.0;
    Timer timer("fixed_subsimW");

    _hyperGraph.BuildWRRsets(num_samples);
    _hyperGraphVldt.BuildWRRsets(num_samples);
    _numRRsets = _hyperGraph.get_RR_sets_size();
    time1 += timer.get_operation_time();

    const auto infSelf = MaxWCover(targetSize);
    time2 += timer.get_operation_time();

    _res.set_running_time(timer.get_total_time());
    _res.set_influence_original(infSelf);
    _res.set_seed_vec(_vecSeed);
    _res.set_RR_sets_size(_numRRsets);
    _res.set_sampling_time(time1);
    _res.set_selection_time(time2);

    std::cout << "==>Time for weighted RR sets and greedy: " << time1 << ", " << time2 << '\n';
    return 0.0;
}

// ===== FRIM xi selection (implement-spec.md) =====

std::vector<double> Alg::buildCumulativeWeights(const std::vector<double>& weights)
{
    std::vector<double> cum(weights.size(), 0.0);
    double sum = 0.0;
    for (size_t i = 0; i < weights.size(); i++)
    {
        sum += std::max(0.0, weights[i]);
        cum[i] = sum;
    }
    return cum;
}

uint32_t Alg::sampleByCumulativeWeights(const std::vector<double>& cumWeights)
{
    if (cumWeights.empty())
        return 0;
    if (cumWeights.back() <= 0.0)
        return dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(cumWeights.size()));

    const double target = dsfmt_gv_genrand_close_open() * cumWeights.back();
    if (target < cumWeights[0])
        return 0;

    int left = 0;
    int right = static_cast<int>(cumWeights.size()) - 1;
    while (left < right)
    {
        const int mid = left + (right - left) / 2;
        if (cumWeights[mid] <= target)
            left = mid + 1;
        else
            right = mid;
    }
    return static_cast<uint32_t>(left);
}

double Alg::nodeObjectiveValue(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi,
    uint32_t nodeId)
{
    return tau[nodeId] * (1.0 - lam[nodeId] * xi[nodeId]);
}

static void frimFillXiStats(FrimXiResult& result, double xi_lo)
{
    size_t num_one = 0;
    for (double x : result.xi)
    {
        if (x > 0.5)
            num_one++;
    }
    result.run_info.num_xi_one = num_one;
    result.run_info.num_xi_lo = result.xi.size() - num_one;
    result.run_info.xi_lo = xi_lo;
}

static std::vector<double> frimRootSamplingWeights(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi)
{
    std::vector<double> weights(tau.size(), 0.0);
    for (size_t v = 0; v < tau.size(); v++)
    {
        if (v < lam.size() && v < xi.size())
            weights[v] = std::max(0.0, tau[v] * (1.0 - lam[v] * xi[v]));
    }
    return weights;
}

static double frimRootWeightSum(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi)
{
    double sum = 0.0;
    for (size_t v = 0; v < tau.size(); v++)
    {
        if (v < lam.size() && v < xi.size())
            sum += std::max(0.0, tau[v] * (1.0 - lam[v] * xi[v]));
    }
    return sum;
}

FrimRRSample Alg::frimSampleOneRR(
    const std::vector<double>& cum_root_weights,
    const std::vector<double>& q,
    const std::vector<double>& xi) const
{
    const uint32_t v_root = sampleByCumulativeWeights(cum_root_weights);
    return _hyperGraph.BuildOneFrimRRSample(v_root, q, xi);
}

FrimRRSample Alg::frimReverseBfsFirstHitRs(
    uint32_t v_root,
    const std::vector<double>& q,
    const std::vector<double>& xi) const
{
    return _hyperGraph.BuildOneFrimRRSample(v_root, q, xi);
}

std::vector<FrimRRSample> Alg::frimBuildRRSamples(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi) const
{
    std::vector<FrimRRSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumRoot =
        buildCumulativeWeights(frimRootSamplingWeights(tau, lam, xi));

    for (size_t r = 0; r < num_rr; r++)
        samples.push_back(frimSampleOneRR(cumRoot, q, xi));
    return samples;
}

double Alg::frimEstimateJ(
    const std::vector<double>& xi,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    size_t num_rr)
{
    const double root_weight_sum = frimRootWeightSum(tau, lam, xi);
    if (num_rr == 0 || root_weight_sum <= 0.0)
        return 0.0;

    _hyperGraph.FrimClearRRTracking();
    const std::vector<double> cumRoot =
        buildCumulativeWeights(frimRootSamplingWeights(tau, lam, xi));
    for (size_t r = 0; r < num_rr; r++)
    {
        const FrimRRSample sample = frimSampleOneRR(cumRoot, q, xi);
        const double weight = sample.hit ? 1.0 : 0.0;
        _hyperGraph.FrimRegisterRRSample(sample.rr_nodes, weight);
    }
    return _hyperGraph.FrimEstimateJFromDummy(root_weight_sum, num_rr);
}

void Alg::frimMethodRR(
    const std::vector<FrimRRSample>& samples,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double root_weight_sum,
    size_t num_rr,
    double xi_lo,
    std::vector<double>& xi,
    std::vector<double>& alpha) const
{
    std::vector<double> K(_numV, 0.0);
    std::vector<double> P(_numV, 0.0);

    for (const auto& sample : samples)
    {
        if (!sample.hit)
            continue;

        const uint32_t v = sample.root;
        const double coef = tau[v] * (1.0 - lam[v]);
        for (uint32_t x : sample.rr_nodes)
        {
            P[x] += 1.0;
            if (x != v)
                K[x] += coef;
        }
    }

    const double scale = root_weight_sum / static_cast<double>(num_rr);
    xi.assign(_numV, xi_lo);
    alpha.assign(_numV, 0.0);
    for (size_t u = 0; u < _numV; u++)
    {
        K[u] *= scale;
        P[u] *= scale;
        alpha[u] = K[u] - tau[u] * lam[u] * P[u];
        xi[u] = (alpha[u] > 0.0) ? 1.0 : xi_lo;
    }
}

void Alg::frimComputeActive(
    const FrimLiveEdgeSample& sample,
    const std::vector<double>& q,
    const std::vector<double>& xi,
    std::vector<uint8_t>& active,
    uint32_t silent_node) const
{
    active.assign(_numV, 0);
    std::queue<uint32_t> que;
    for (size_t u = 0; u < _numV; u++)
    {
        if (sample.seed_gates[u] > q[u])
            continue;
        active[u] = 1;
        que.push(static_cast<uint32_t>(u));
    }

    while (!que.empty())
    {
        const uint32_t u = que.front();
        que.pop();

        if (u == silent_node)
            continue;
        if (sample.gates[u] > xi[u])
            continue;

        for (uint32_t v : sample.kept_out[u])
        {
            if (active[v])
                continue;
            active[v] = 1;
            que.push(v);
        }
    }
}

std::vector<uint8_t> Alg::frimForwardReachFrom(
    const FrimLiveEdgeSample& sample,
    const std::vector<double>& xi,
    uint32_t src) const
{
    std::vector<uint8_t> reached(_numV, 0);
    if (src >= _numV || sample.gates[src] > xi[src])
        return reached;

    std::queue<uint32_t> q;
    reached[src] = 1;
    q.push(src);
    while (!q.empty())
    {
        const uint32_t u = q.front();
        q.pop();
        if (sample.gates[u] > xi[u])
            continue;

        for (uint32_t v : sample.kept_out[u])
        {
            if (reached[v])
                continue;
            reached[v] = 1;
            q.push(v);
        }
    }
    return reached;
}

double Alg::frimObjectiveOnSample(
    const FrimLiveEdgeSample& sample,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi) const
{
    double sum = 0.0;
    for (size_t v = 0; v < _numV; v++)
    {
        if (sample.active[v])
            sum += nodeObjectiveValue(tau, lam, xi, static_cast<uint32_t>(v));
    }
    return sum;
}

double Alg::frimFlipDelta(
    const FrimLiveEdgeSample& sample,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi,
    uint32_t u,
    double xi_new) const
{
    if (u >= _numV || !sample.active[u])
        return 0.0;

    const double xi_old = xi[u];
    if (std::abs(xi_old - xi_new) <= 0.0)
        return 0.0;

    double delta = 0.0;
    if (xi_old > xi_new)
    {
        std::vector<uint8_t> active_minus;
        frimComputeActive(sample, q, xi, active_minus, u);

        for (size_t v = 0; v < _numV; v++)
        {
            if (sample.active[v] && !active_minus[v])
                delta -= nodeObjectiveValue(tau, lam, xi, static_cast<uint32_t>(v));
        }
        delta += tau[u] * lam[u] * (xi_old - xi_new);
    }
    else
    {
        std::vector<double> xi_eff = xi;
        xi_eff[u] = xi_new;
        const std::vector<uint8_t> reached = frimForwardReachFrom(sample, xi_eff, u);
        for (size_t v = 0; v < _numV; v++)
        {
            if (reached[v] && !sample.active[v])
                delta += nodeObjectiveValue(tau, lam, xi, static_cast<uint32_t>(v));
        }
        delta += tau[u] * lam[u] * (xi_old - xi_new);
    }
    return delta;
}

void Alg::frimApplyFlip(
    FrimLiveEdgeSample& sample,
    const std::vector<double>& q,
    const std::vector<double>& xi,
    uint32_t u,
    double xi_new,
    double xi_lo) const
{
    if (u >= _numV || !sample.active[u])
        return;

    if (std::abs(xi_new - xi_lo) <= 0.0)
    {
        frimComputeActive(sample, q, xi, sample.active, u);
    }
    else
    {
        const std::vector<uint8_t> reached = frimForwardReachFrom(sample, xi, u);
        for (size_t v = 0; v < _numV; v++)
        {
            if (reached[v])
                sample.active[v] = 1;
        }
    }
}

FrimXiResult Alg::frimMethodMC(
    const Graph& forwardGraph,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_mc,
    int max_sweeps,
    double delta_tol) const
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_mc == 0 || q.size() != _numV)
        return result;

    Timer timer;
    std::vector<FrimLiveEdgeSample> live_samples;
    live_samples.reserve(num_mc);
    for (size_t r = 0; r < num_mc; r++)
    {
        FrimLiveEdgeSample sample;
        sample.kept_out.assign(forwardGraph.size(), {});
        sample.gates.assign(_numV, 0.0);
        for (size_t u = 0; u < forwardGraph.size(); u++)
        {
            for (const auto& edge : forwardGraph[u])
            {
                if (dsfmt_gv_genrand_open_close() <= edge.second)
                    sample.kept_out[u].push_back(edge.first);
            }
        }
        for (size_t i = 0; i < _numV; i++)
            sample.gates[i] = dsfmt_gv_genrand_close_open();
        sample.seed_gates.assign(_numV, 0.0);
        for (size_t i = 0; i < _numV; i++)
            sample.seed_gates[i] = dsfmt_gv_genrand_open_close();

        frimComputeActive(sample, q, result.xi, sample.active);
        live_samples.push_back(std::move(sample));
    }
    result.run_info.time_sample_sec = timer.get_operation_time();

    double J_hat = 0.0;
    for (const auto& sample : live_samples)
        J_hat += frimObjectiveOnSample(sample, tau, lam, result.xi);
    J_hat /= static_cast<double>(num_mc);

    std::vector<uint32_t> order(_numV);
    std::iota(order.begin(), order.end(), 0U);
    int sweeps_completed = 0;
    for (int sweep = 0; sweep < max_sweeps; sweep++)
    {
        sweeps_completed++;
        bool improved = false;
        for (size_t i = _numV - 1; i > 0; i--)
        {
            const size_t j = dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
            std::swap(order[i], order[j]);
        }

        for (uint32_t u : order)
        {
            const double xi_new = (result.xi[u] > 0.5) ? xi_lo : 1.0;
            if (std::abs(result.xi[u] - xi_new) <= 0.0)
                continue;

            double delta = 0.0;
            for (const auto& sample : live_samples)
                delta += frimFlipDelta(sample, q, tau, lam, result.xi, u, xi_new);
            delta /= static_cast<double>(num_mc);

            if (delta > delta_tol)
            {
                result.xi[u] = xi_new;
                for (auto& sample : live_samples)
                    frimApplyFlip(sample, q, result.xi, u, xi_new, xi_lo);
                J_hat += delta;
                improved = true;
            }
        }

        if (!improved)
            break;
    }
    result.run_info.time_solve_sec = timer.get_operation_time();
    result.run_info.sweeps_completed = sweeps_completed;
    result.run_info.max_sweeps = max_sweeps;

    result.J_method_mc = J_hat;
    result.run_info.time_total_sec = timer.get_total_time();
    return result;
}

FrimXiResult Alg::frim_solve_rr(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    double tau_sum = 0.0;
    for (double t : tau)
        tau_sum += t;
    if (tau_sum <= 0.0)
        return result;

    std::cout << "[FRIM-RR] sampling " << num_rr
              << " tau*(1-lambda*xi)-weighted RR paths (random seed)..." << std::endl;
    Timer timer;
    const std::vector<double> xi_all_one(_numV, 1.0);
    const std::vector<FrimRRSample> rr_samples =
        frimBuildRRSamples(num_rr, tau, lam, q, xi_all_one);
    const double root_weight_sum = frimRootWeightSum(tau, lam, xi_all_one);
    result.run_info.time_sample_sec = timer.get_operation_time();
    frimMethodRR(rr_samples, tau, lam, q, root_weight_sum, num_rr, xi_lo, result.xi, result.alpha);
    result.run_info.time_solve_sec = timer.get_operation_time();
    result.J_method_rr = frimEstimateJ(result.xi, q, tau, lam, num_rr);
    result.run_info.time_estimate_sec = timer.get_operation_time();
    result.J_hat = result.J_method_rr;
    result.run_info.time_total_sec = timer.get_total_time();
    std::cout << "[FRIM-RR] done: J_hat=" << result.J_method_rr << std::endl;

    frimFillXiStats(result, xi_lo);
    size_t num_xi_one = result.run_info.num_xi_one;
    std::cout << "[FRIM-RR] xi: " << num_xi_one << " nodes at 1, "
              << (_numV - num_xi_one) << " nodes at " << xi_lo << std::endl;
    return result;
}

FrimXiResult Alg::frimMethodRRNaive(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_rr == 0 || q.size() != _numV
        || tau.size() != _numV || lam.size() != _numV)
        return result;

    const double root_weight_sum = frimRootWeightSum(tau, lam, result.xi);
    if (root_weight_sum <= 0.0)
        return result;

    Timer timer;
    // Estimator: tau*(1-lam*xi) root sampling, hit weight 1; fresh R samples per J call.
    result.J_hat = frimEstimateJ(result.xi, q, tau, lam, num_rr);
    result.run_info.time_sample_sec = timer.get_operation_time();
    std::cout << "[FRIM-RR-naive] initial J_hat=" << result.J_hat
              << " (tau*(1-lam*xi) root, hit=1, R=" << num_rr << ")"
              << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

    int sweeps_completed = 0;
    for (int sweep = 0; sweep < max_sweeps; sweep++)
    {
        sweeps_completed++;
        bool improved = false;
        std::vector<uint32_t> order(_numV);
        std::iota(order.begin(), order.end(), 0U);
        for (size_t i = _numV - 1; i > 0; i--)
        {
            const size_t j = dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
            std::swap(order[i], order[j]);
        }

        std::cout << "[FRIM-RR-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " start, J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        uint32_t currProgress = 0;
        size_t node_idx = 0;
        size_t updates_this_sweep = 0;
        for (uint32_t u : order)
        {
            node_idx++;
            if (node_idx * 100 >= _numV * currProgress)
            {
                std::cout << "[FRIM-RR-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                          << " nodes " << currProgress << "%, updates=" << updates_this_sweep
                          << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
                currProgress += 20;
            }

            std::vector<double> xi_at_lo = result.xi;
            std::vector<double> xi_at_one = result.xi;
            xi_at_lo[u] = xi_lo;
            xi_at_one[u] = 1.0;

            const double J_lo = frimEstimateJ(xi_at_lo, q, tau, lam, num_rr);
            const double J_one = frimEstimateJ(xi_at_one, q, tau, lam, num_rr);

            const double xi_pick =
                (J_one > J_lo + delta_tol) ? 1.0
                : (J_lo > J_one + delta_tol) ? xi_lo
                : (J_one >= J_lo ? 1.0 : xi_lo);

            if (std::abs(result.xi[u] - xi_pick) > 0.0)
            {
                result.xi[u] = xi_pick;
                improved = true;
                updates_this_sweep++;
            }
        }

        result.J_hat = frimEstimateJ(result.xi, q, tau, lam, num_rr);
        std::cout << "[FRIM-RR-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " done: updates=" << updates_this_sweep
                  << ", improved=" << (improved ? "yes" : "no")
                  << ", J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        if (!improved)
            break;
    }

    result.run_info.time_solve_sec = timer.get_operation_time();
    result.run_info.sweeps_completed = sweeps_completed;
    result.run_info.max_sweeps = max_sweeps;
    result.J_method_rr_naive = frimEstimateJ(result.xi, q, tau, lam, num_rr);
    result.J_hat = result.J_method_rr_naive;
    result.run_info.time_total_sec = timer.get_total_time();
    return result;
}

static double frimTauSum(const std::vector<double>& tau)
{
    double sum = 0.0;
    for (double t : tau)
        sum += std::max(0.0, t);
    return sum;
}

FrimRRStructureSample Alg::frimBuildOneRRStructure(
    const std::vector<double>& cum_tau,
    const std::vector<double>& q) const
{
    const uint32_t v_root = sampleByCumulativeWeights(cum_tau);
    return _hyperGraph.BuildOneFrimRRStructure(v_root, q);
}

std::vector<FrimRRStructureSample> Alg::frimBuildRRStructureSamples(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q) const
{
    std::vector<FrimRRStructureSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    for (size_t r = 0; r < num_rr; r++)
        samples.push_back(frimBuildOneRRStructure(cumTau, q));
    return samples;
}

double Alg::frimRRStructureSampleWeight(
    const FrimRRStructureSample& sample,
    const std::vector<double>& xi,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    uint32_t override_node,
    double override_xi)
{
    if (!sample.hit)
        return 0.0;

    const auto xiAt = [&](uint32_t v) -> double
    {
        if (v == override_node)
            return override_xi;
        return (v < xi.size()) ? xi[v] : 0.0;
    };

    const size_t numV = sample.parent.size();
    std::vector<uint8_t> reachable(numV, 0);

    for (uint32_t u : sample.bfs_order)
    {
        if (u == sample.root)
        {
            reachable[u] = 1;
        }
        else
        {
            const uint32_t p = sample.parent[u];
            if (p == UINT32_MAX || !reachable[p])
                continue;
            if (sample.xi_gate[u] < 0.0f)
                continue;
            if (u != sample.root && static_cast<double>(sample.xi_gate[u]) > xiAt(u))
                continue;
            reachable[u] = 1;
        }

        if (reachable[u] && u == sample.hit_node)
        {
            double log_w = 0.0;
            uint32_t cur = sample.hit_node;
            while (cur != sample.root)
            {
                if (cur != sample.hit_node && cur != sample.root)
                {
                    const double xv = xiAt(cur);
                    if (xv <= 0.0)
                        return 0.0;
                    log_w += std::log(xv);
                }
                const uint32_t p = sample.parent[cur];
                if (p == UINT32_MAX)
                    return 0.0;
                cur = p;
            }

            if (sample.root >= tau.size() || sample.root >= lam.size())
                return 0.0;

            const double root_coef = 1.0 - lam[sample.root] * xiAt(sample.root);
            return root_coef * std::exp(log_w);
        }
    }

    return 0.0;
}

double Alg::frimEstimateJFromRRStructures(
    const std::vector<FrimRRStructureSample>& samples,
    const std::vector<double>& xi,
    const std::vector<double>& tau,
    const std::vector<double>& lam) const
{
    if (samples.empty())
        return 0.0;

    const double tau_sum = frimTauSum(tau);
    if (tau_sum <= 0.0)
        return 0.0;

    double total = 0.0;
    for (const auto& sample : samples)
        total += frimRRStructureSampleWeight(sample, xi, tau, lam);

    return (tau_sum / static_cast<double>(samples.size())) * total;
}

FrimXiResult Alg::frimMethodRRCrn(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_rr == 0 || q.size() != _numV
        || tau.size() != _numV || lam.size() != _numV)
        return result;

    const double tau_sum = frimTauSum(tau);
    if (tau_sum <= 0.0)
        return result;

    Timer timer;
    const std::vector<FrimRRStructureSample> samples =
        frimBuildRRStructureSamples(num_rr, tau, q);
    result.run_info.time_sample_sec = timer.get_operation_time();

    std::vector<std::vector<size_t>> node_samples(_numV);
    for (size_t r = 0; r < samples.size(); r++)
    {
        for (uint32_t u : samples[r].bfs_order)
            node_samples[u].push_back(r);
    }

    std::vector<double> sample_weights(num_rr, 0.0);
    double weight_sum = 0.0;
    for (size_t r = 0; r < samples.size(); r++)
    {
        sample_weights[r] =
            frimRRStructureSampleWeight(samples[r], result.xi, tau, lam);
        weight_sum += sample_weights[r];
    }

    const double scale = tau_sum / static_cast<double>(num_rr);
    result.J_hat = scale * weight_sum;
    std::cout << "[FRIM-RR-CRN] initial J_hat=" << result.J_hat
              << " (tau-root, xi_gate CRN, R=" << num_rr << ")"
              << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

    auto recomputeAffected = [&](uint32_t u, double xi_u, double& affected_sum)
    {
        affected_sum = 0.0;
        for (size_t r : node_samples[u])
        {
            affected_sum += frimRRStructureSampleWeight(
                samples[r], result.xi, tau, lam, u, xi_u);
        }
    };

    int sweeps_completed = 0;
    for (int sweep = 0; sweep < max_sweeps; sweep++)
    {
        sweeps_completed++;
        bool improved = false;
        std::vector<uint32_t> order(_numV);
        std::iota(order.begin(), order.end(), 0U);
        for (size_t i = _numV - 1; i > 0; i--)
        {
            const size_t j = dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
            std::swap(order[i], order[j]);
        }

        std::cout << "[FRIM-RR-CRN] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " start, J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        uint32_t currProgress = 0;
        size_t node_idx = 0;
        size_t updates_this_sweep = 0;
        for (uint32_t u : order)
        {
            node_idx++;
            if (node_idx * 100 >= _numV * currProgress)
            {
                std::cout << "[FRIM-RR-CRN] sweep " << (sweep + 1) << "/" << max_sweeps
                          << " nodes " << currProgress << "%, updates=" << updates_this_sweep
                          << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
                currProgress += 20;
            }

            double affected_old = 0.0;
            for (size_t r : node_samples[u])
                affected_old += sample_weights[r];

            double affected_lo = 0.0;
            double affected_one = 0.0;
            recomputeAffected(u, xi_lo, affected_lo);
            recomputeAffected(u, 1.0, affected_one);

            const double J_lo = scale * (weight_sum - affected_old + affected_lo);
            const double J_one = scale * (weight_sum - affected_old + affected_one);

            const double xi_pick =
                (J_one > J_lo + delta_tol) ? 1.0
                : (J_lo > J_one + delta_tol) ? xi_lo
                : (J_one >= J_lo ? 1.0 : xi_lo);

            if (std::abs(result.xi[u] - xi_pick) > 0.0)
            {
                result.xi[u] = xi_pick;
                const double affected_new =
                    (std::abs(xi_pick - 1.0) <= 0.0) ? affected_one : affected_lo;
                weight_sum = weight_sum - affected_old + affected_new;
                for (size_t r : node_samples[u])
                {
                    sample_weights[r] = frimRRStructureSampleWeight(
                        samples[r], result.xi, tau, lam);
                }
                result.J_hat = scale * weight_sum;
                improved = true;
                updates_this_sweep++;
            }
        }

        std::cout << "[FRIM-RR-CRN] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " done: updates=" << updates_this_sweep
                  << ", improved=" << (improved ? "yes" : "no")
                  << ", J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        if (!improved)
            break;
    }

    result.run_info.time_solve_sec = timer.get_operation_time();
    result.run_info.sweeps_completed = sweeps_completed;
    result.run_info.max_sweeps = max_sweeps;
    result.J_method_rr_crn = scale * weight_sum;
    result.J_hat = result.J_method_rr_crn;
    result.run_info.time_total_sec = timer.get_total_time();
    return result;
}

FrimXiResult Alg::frim_solve_rr_crn(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-RR-CRN] " << num_rr
              << " tau-root structural RR samples, CRN {xi_lo,1} per-node..."
              << std::endl;
    const FrimXiResult method_crn = frimMethodRRCrn(
        tau, lam, q, init, xi_lo, num_rr, max_sweeps, delta_tol);
    result = method_crn;
    result.J_hat = result.J_method_rr_crn;
    std::cout << "[FRIM-RR-CRN] done: J_hat=" << result.J_method_rr_crn << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-RR-CRN] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

double Alg::frim_estimate_j_at_xi(
    const std::vector<double>& xi,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    size_t num_rr)
{
    return frimEstimateJ(xi, q, tau, lam, num_rr);
}

double Alg::frim_estimate_mc_gate_at_xi(
    const Graph& forwardGraph,
    const std::vector<double>& xi,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    size_t num_mc) const
{
    if (num_mc == 0)
        return 0.0;

    double sum = 0.0;
    for (size_t r = 0; r < num_mc; r++)
    {
        FrimLiveEdgeSample sample = frimSampleLiveEdge(forwardGraph);
        frimComputeActive(sample, q, xi, sample.active);
        sum += frimObjectiveOnSample(sample, tau, lam, xi);
    }
    return sum / static_cast<double>(num_mc);
}

FrimXiResult Alg::frim_solve_rr_naive(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-RR-naive] " << num_rr
              << " tau*(1-lam*xi)-root RR samples, per-node {xi_lo,1} J estimate, resample per call..."
              << std::endl;
    const FrimXiResult method_naive = frimMethodRRNaive(
        tau, lam, q, init, xi_lo, num_rr, max_sweeps, delta_tol);
    result = method_naive;
    result.J_hat = result.J_method_rr_naive;
    std::cout << "[FRIM-RR-naive] done: J_hat=" << result.J_method_rr_naive << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-RR-naive] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

FrimLiveEdgeSample Alg::frimSampleLiveEdge(
    const Graph& forwardGraph) const
{
    FrimLiveEdgeSample sample;
    sample.kept_out.assign(forwardGraph.size(), {});
    sample.gates.assign(_numV, 0.0);
    sample.seed_gates.assign(_numV, 0.0);
    for (size_t u = 0; u < forwardGraph.size(); u++)
    {
        for (const auto& edge : forwardGraph[u])
        {
            if (dsfmt_gv_genrand_open_close() <= edge.second)
                sample.kept_out[u].push_back(edge.first);
        }
    }
    for (size_t i = 0; i < _numV; i++)
    {
        sample.gates[i] = dsfmt_gv_genrand_close_open();
        sample.seed_gates[i] = dsfmt_gv_genrand_open_close();
    }
    return sample;
}

double Alg::frimObjectiveWithXi(
    const FrimLiveEdgeSample& sample,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& xi) const
{
    std::vector<uint8_t> active;
    frimComputeActive(sample, q, xi, active);

    double sum = 0.0;
    for (size_t v = 0; v < _numV; v++)
    {
        if (active[v])
            sum += nodeObjectiveValue(tau, lam, xi, static_cast<uint32_t>(v));
    }
    return sum;
}

FrimXiResult Alg::frimMethodMCNaive(
    const Graph& forwardGraph,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_mc,
    int max_sweeps,
    double delta_tol) const
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_mc == 0 || q.size() != _numV)
        return result;

    Timer timer;
    uint32_t currProgress = 0;

    std::vector<FrimLiveEdgeSample> live_samples;
    live_samples.reserve(num_mc);
    for (size_t r = 0; r < num_mc; r++)
    {
        if (r * 100 >= num_mc * currProgress)
        {
            std::cout << "[FRIM-MC-naive] sampling live-edge: " << currProgress << "%, "
                      << "elapsed=" << timer.get_total_time() << "s" << std::endl;
            currProgress += 20;
        }
        live_samples.push_back(frimSampleLiveEdge(forwardGraph));
    }
    std::cout << "[FRIM-MC-naive] sampled " << num_mc << " live-edge graphs, "
              << "elapsed=" << timer.get_total_time() << "s" << std::endl;
    result.run_info.time_sample_sec = timer.get_operation_time();

    int sweeps_completed = 0;
    for (int sweep = 0; sweep < max_sweeps; sweep++)
    {
        sweeps_completed++;
        bool improved = false;
        size_t flips_this_sweep = 0;
        std::vector<uint32_t> order(_numV);
        std::iota(order.begin(), order.end(), 0U);
        for (size_t i = _numV - 1; i > 0; i--)
        {
            const size_t j = dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
            std::swap(order[i], order[j]);
        }

        std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " start, elapsed=" << timer.get_total_time() << "s" << std::endl;

        currProgress = 0;
        size_t node_idx = 0;
        for (uint32_t u : order)
        {
            node_idx++;
            if (node_idx * 100 >= _numV * currProgress)
            {
                std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                          << " nodes " << currProgress << "%, flips=" << flips_this_sweep
                          << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
                currProgress += 20;
            }

            const double xi_new = (result.xi[u] > 0.5) ? xi_lo : 1.0;
            if (std::abs(result.xi[u] - xi_new) <= 0.0)
                continue;

            std::vector<double> xi_flip = result.xi;
            xi_flip[u] = xi_new;

            double sum_orig = 0.0;
            double sum_flip = 0.0;
            for (const auto& sample : live_samples)
            {
                sum_orig += frimObjectiveWithXi(sample, q, tau, lam, result.xi);
                sum_flip += frimObjectiveWithXi(sample, q, tau, lam, xi_flip);
            }
            const double delta = (sum_flip - sum_orig) / static_cast<double>(num_mc);

            if (delta > delta_tol)
            {
                result.xi[u] = xi_new;
                improved = true;
                flips_this_sweep++;
            }
        }

        std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                  << " done: flips=" << flips_this_sweep
                  << ", improved=" << (improved ? "yes" : "no")
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        if (!improved)
            break;
    }
    result.run_info.time_solve_sec = timer.get_operation_time();
    result.run_info.sweeps_completed = sweeps_completed;
    result.run_info.max_sweeps = max_sweeps;

    std::cout << "[FRIM-MC-naive] evaluating J_hat..." << std::endl;
    double J_hat = 0.0;
    for (const auto& sample : live_samples)
        J_hat += frimObjectiveWithXi(sample, q, tau, lam, result.xi);
    result.J_method_mc_naive = J_hat / static_cast<double>(num_mc);
    result.J_hat = result.J_method_mc_naive;
    result.run_info.time_estimate_sec = timer.get_operation_time();
    result.run_info.time_total_sec = timer.get_total_time();
    std::cout << "[FRIM-MC-naive] J_hat evaluated, elapsed=" << result.run_info.time_total_sec << "s"
              << std::endl;
    return result;
}

FrimXiResult Alg::frim_solve_mc(
    const Graph& forwardGraph,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_mc,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-MC] coordinate ascent with " << num_mc << " CRN samples (random seed)..." << std::endl;
    const FrimXiResult method_mc = frimMethodMC(
        forwardGraph, tau, lam, q, init, xi_lo, num_mc, max_sweeps, delta_tol);
    result = method_mc;
    result.J_hat = result.J_method_mc;
    std::cout << "[FRIM-MC] done: J_hat=" << result.J_method_mc << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-MC] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

FrimXiResult Alg::frim_solve_mc_naive(
    const Graph& forwardGraph,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_mc,
    int max_sweeps,
    double delta_tol)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-MC-naive] " << num_mc
              << " CRN live-edge samples (full BFS delta, real-time xi)..." << std::endl;
    const FrimXiResult method_naive = frimMethodMCNaive(
        forwardGraph, tau, lam, q, init, xi_lo, num_mc, max_sweeps, delta_tol);
    result = method_naive;
    result.J_hat = result.J_method_mc_naive;
    std::cout << "[FRIM-MC-naive] done: J_hat=" << result.J_method_mc_naive << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-MC-naive] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

// ===== end FRIM xi selection =====
