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

static double frimTauSum(const std::vector<double>& tau)
{
    double sum = 0.0;
    for (double t : tau)
        sum += std::max(0.0, t);
    return sum;
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

std::vector<FrimRRSample> Alg::frimBuildTauRootRRSamplesAtXiOne(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q) const
{
    std::vector<FrimRRSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    const std::vector<double> xi_one(_numV, 1.0);
    for (size_t r = 0; r < num_rr; r++)
    {
        const uint32_t v_root = sampleByCumulativeWeights(cumTau);
        samples.push_back(_hyperGraph.BuildOneFrimRRSample(v_root, q, xi_one));
    }
    return samples;
}

std::vector<FrimRRGraphSample> Alg::frimBuildTauRootRRGraphSamplesAtXiOne(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q) const
{
    std::vector<FrimRRGraphSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    for (size_t r = 0; r < num_rr; r++)
    {
        const uint32_t v_root = sampleByCumulativeWeights(cumTau);
        samples.push_back(_hyperGraph.BuildOneFrimRRGraphSample(v_root, q));
    }
    return samples;
}

std::vector<FrimRRGraphSample> Alg::frimBuildTauRootRRGraphSamples(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q,
    const std::vector<double>& xi) const
{
    std::vector<FrimRRGraphSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    for (size_t r = 0; r < num_rr; r++)
    {
        const uint32_t v_root = sampleByCumulativeWeights(cumTau);
        samples.push_back(_hyperGraph.BuildOneFrimRRGraphSampleWithXi(v_root, q, xi));
    }
    return samples;
}

std::vector<FrimRRGraphSample> Alg::frimBuildTauRootRRGraphSamplesCrn(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q,
    bool store_hit_only,
    size_t* num_discarded_out) const
{
    std::vector<FrimRRGraphSample> samples;
    samples.reserve(num_rr);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    size_t num_discarded = 0;
    for (size_t r = 0; r < num_rr; r++)
    {
        const uint32_t v_root = sampleByCumulativeWeights(cumTau);
        FrimRRGraphSample sample =
            _hyperGraph.BuildOneFrimRRGraphSampleCrn(v_root, q);
        if (store_hit_only && !sample.hit)
        {
            // Keep a lightweight placeholder so samples.size() == num_rr (scale T/R).
            num_discarded++;
            FrimRRGraphSample ph;
            ph.root = v_root;
            ph.hit = false;
            samples.push_back(std::move(ph));
            continue;
        }
        samples.push_back(std::move(sample));
    }
    if (num_discarded_out != nullptr)
        *num_discarded_out = num_discarded;
    return samples;
}

namespace
{
struct FrimRootCountSummary
{
    size_t min_count = 0;
    size_t max_count = 0;
    double mean_count = 0.0;
    size_t num_at_min = 0;
    std::vector<uint32_t> min_nodes;
};

static FrimRootCountSummary frimSummarizeRootCounts(
    const std::vector<size_t>& root_counts,
    size_t max_examples = 5)
{
    FrimRootCountSummary summary;
    if (root_counts.empty())
        return summary;

    size_t total = 0;
    summary.min_count = root_counts[0];
    summary.max_count = root_counts[0];
    for (size_t u = 0; u < root_counts.size(); u++)
    {
        const size_t c = root_counts[u];
        total += c;
        summary.min_count = std::min(summary.min_count, c);
        summary.max_count = std::max(summary.max_count, c);
    }
    summary.mean_count =
        static_cast<double>(total) / static_cast<double>(root_counts.size());

    for (size_t u = 0; u < root_counts.size(); u++)
    {
        if (root_counts[u] != summary.min_count)
            continue;
        summary.num_at_min++;
        if (summary.min_nodes.size() < max_examples)
            summary.min_nodes.push_back(static_cast<uint32_t>(u));
    }
    return summary;
}

static void frimPrintRootCountSummary(
    const char* label,
    const FrimRootCountSummary& summary)
{
    std::cout << "[FRIM-RR-ROOT-STAT] " << label
              << ": min=" << summary.min_count
              << ", max=" << summary.max_count
              << ", mean=" << summary.mean_count
              << ", num_nodes_at_min=" << summary.num_at_min;
    if (!summary.min_nodes.empty())
    {
        std::cout << ", example_min_nodes={";
        for (size_t i = 0; i < summary.min_nodes.size(); i++)
        {
            if (i > 0)
                std::cout << ',';
            std::cout << summary.min_nodes[i];
        }
        std::cout << '}';
    }
    std::cout << std::endl;
}
}  // namespace

void Alg::frimStatAndPrintRRGraphRootCounts(
    size_t num_rr,
    const std::vector<double>& tau,
    const std::vector<double>& q,
    bool store_hit_only,
    bool build_rr_graphs) const
{
    if (_numV == 0 || num_rr == 0 || tau.size() != _numV || q.size() != _numV)
    {
        std::cout << "[FRIM-RR-ROOT-STAT] invalid inputs: num_v=" << _numV
                  << " num_rr=" << num_rr << " tau=" << tau.size()
                  << " q=" << q.size() << std::endl;
        return;
    }

    const double tau_sum = frimTauSum(tau);
    std::cout << "[FRIM-RR-ROOT-STAT] R=" << num_rr
              << " num_v=" << _numV
              << " tau_sum=" << tau_sum
              << " store_hit_only=" << (store_hit_only ? "on" : "off")
              << " build_rr_graphs=" << (build_rr_graphs ? "on" : "off")
              << std::endl;

    std::vector<size_t> root_draws(_numV, 0);
    std::vector<size_t> root_stored(_numV, 0);
    const std::vector<double> cumTau = buildCumulativeWeights(tau);
    size_t num_discarded = 0;
    size_t num_stored = 0;

    Timer timer;
    for (size_t r = 0; r < num_rr; r++)
    {
        const uint32_t v_root = sampleByCumulativeWeights(cumTau);
        root_draws[v_root]++;

        if (!build_rr_graphs)
            continue;

        const FrimRRGraphSample sample =
            _hyperGraph.BuildOneFrimRRGraphSampleCrn(v_root, q);
        if (store_hit_only && !sample.hit)
        {
            num_discarded++;
            continue;
        }
        root_stored[v_root]++;
        num_stored++;
    }

    const FrimRootCountSummary draw_summary =
        frimSummarizeRootCounts(root_draws);
    frimPrintRootCountSummary("tau-root draws (all R attempts)", draw_summary);

    if (build_rr_graphs)
    {
        const FrimRootCountSummary stored_summary =
            frimSummarizeRootCounts(root_stored);
        std::cout << "[FRIM-RR-ROOT-STAT] stored_rr_graphs=" << num_stored
                  << "/" << num_rr;
        if (store_hit_only)
            std::cout << ", discarded_no_hit=" << num_discarded;
        std::cout << std::endl;
        frimPrintRootCountSummary("stored RR-graph roots", stored_summary);
    }

    std::cout << "[FRIM-RR-ROOT-STAT] elapsed="
              << timer.get_operation_time() << "s" << std::endl;
}


namespace
{
inline bool frimByteMapGet(const FrimByteMap& m, uint32_t u)
{
    const auto it = m.find(u);
    return it != m.end() && it->second != 0;
}

inline int frimIntMapGet(const FrimIntMap& m, uint32_t u, int def_val = 0)
{
    const auto it = m.find(u);
    return it == m.end() ? def_val : it->second;
}

/// BFS on sample.forward_adj; edge u->v is used only when node_active[v] is set.
void frimFillReachFromSparse(
    uint32_t src,
    const FrimSparseFwdAdj& adj,
    const FrimByteMap& node_active,
    FrimByteMap& reach)
{
    if (!frimByteMapGet(node_active, src))
        return;

    std::queue<uint32_t> que;
    reach[src] = 1;
    que.push(src);

    while (!que.empty())
    {
        const uint32_t u = que.front();
        que.pop();

        const auto it = adj.find(u);
        if (it == adj.end())
            continue;

        for (uint32_t v : it->second)
        {
            if (frimByteMapGet(reach, v) || !frimByteMapGet(node_active, v))
                continue;
            reach[v] = 1;
            que.push(v);
        }
    }
}

void frimFillReachFromMultiSource(
    const std::vector<uint32_t>& sources,
    const FrimSparseFwdAdj& adj,
    const FrimByteMap& node_active,
    FrimByteMap& reach)
{
    reach.clear();
    std::queue<uint32_t> que;
    for (uint32_t src : sources)
    {
        if (!frimByteMapGet(node_active, src))
            continue;
        if (!frimByteMapGet(reach, src))
        {
            reach[src] = 1;
            que.push(src);
        }
    }

    while (!que.empty())
    {
        const uint32_t u = que.front();
        que.pop();

        const auto it = adj.find(u);
        if (it == adj.end())
            continue;

        for (uint32_t v : it->second)
        {
            if (frimByteMapGet(reach, v) || !frimByteMapGet(node_active, v))
                continue;
            reach[v] = 1;
            que.push(v);
        }
    }
}

void frimBuildUndirSparseOnReachable(
    const FrimSparseFwdAdj& forward_adj,
    const FrimByteMap& from_reach,
    const FrimByteMap& node_active,
    FrimSparseFwdAdj& undir_out)
{
    undir_out.clear();
    for (const auto& kv : forward_adj)
    {
        const uint32_t a = kv.first;
        if (!frimByteMapGet(from_reach, a) || !frimByteMapGet(node_active, a))
            continue;

        for (uint32_t b : kv.second)
        {
            if (!frimByteMapGet(from_reach, b) || !frimByteMapGet(node_active, b))
                continue;
            undir_out[a].push_back(b);
            undir_out[b].push_back(a);
        }
    }
}

void frimPlatformTarjanDfsSparse(
    const uint32_t u,
    const FrimSparseFwdAdj& adj,
    FrimIntMap& disc,
    FrimIntMap& low,
    FrimIntMap& parent,
    FrimByteMap& subtree_has_root,
    FrimByteMap& blocks,
    uint32_t root_id,
    int& timer)
{
    disc[u] = low[u] = ++timer;
    subtree_has_root[u] = (u == root_id) ? 1 : 0;

    int children = 0;
    const auto it = adj.find(u);
    if (it == adj.end())
        return;

    const int parent_u = frimIntMapGet(parent, u, -1);
    for (const uint32_t v : it->second)
    {
        if (frimIntMapGet(disc, v) == 0)
        {
            children++;
            parent[v] = static_cast<int>(u);
            frimPlatformTarjanDfsSparse(
                v, adj, disc, low, parent, subtree_has_root, blocks,
                root_id, timer);
            subtree_has_root[u] =
                static_cast<uint8_t>(subtree_has_root[u] | subtree_has_root[v]);
            low[u] = std::min(low[u], low[v]);

            if (parent_u != -1
                && low[v] >= disc[u]
                && subtree_has_root[v])
            {
                blocks[u] = 1;
            }
        }
        else if (static_cast<int>(v) != parent_u)
        {
            low[u] = std::min(low[u], disc[v]);
        }
    }

    if (parent_u == -1 && children > 1)
    {
        for (const uint32_t v : it->second)
        {
            if (frimIntMapGet(parent, v, -1) != static_cast<int>(u))
                continue;
            if (low[v] >= disc[u] && subtree_has_root[v])
                blocks[u] = 1;
        }
    }
}

    // Thread-local sparse buffers (only touched nodes stored)
    thread_local FrimByteMap g_xi_active;
    thread_local FrimByteMap g_from_hits;
    thread_local FrimSparseFwdAdj g_undir;
    thread_local FrimIntMap g_disc;
    thread_local FrimIntMap g_low;
    thread_local FrimIntMap g_parent;
    thread_local FrimByteMap g_subtree_has_root;
    thread_local FrimByteMap g_blocks_tarjan;
    thread_local FrimByteMap g_blocks_xi1;
    thread_local FrimByteMap g_xi_active_xilo;
    thread_local FrimByteMap g_from_hits_xilo;

static const std::vector<uint32_t>& frimRRGraphHitList(
    const FrimRRGraphSample& sample)
{
    return sample.hit_nodes;
}

static bool frimRRGraphSampleHasHits(const FrimRRGraphSample& sample)
{
    return sample.hit && !sample.hit_nodes.empty();
}

bool frimRRGraphRunBlockPipeline(
    const FrimRRGraphSample& sample,
    FrimByteMap& blocks,
    FrimByteMap& xi_active)
{
    blocks.clear();

    const std::vector<uint32_t>& hits = frimRRGraphHitList(sample);
    if (hits.empty() || sample.platform_id == UINT32_MAX)
        return false;

    const uint32_t virtual_seed = sample.platform_id;

    g_from_hits.clear();
    frimFillReachFromMultiSource(
        hits, sample.forward_adj, xi_active, g_from_hits);
    if (!frimByteMapGet(g_from_hits, sample.root))
        return false;

    frimBuildUndirSparseOnReachable(
        sample.forward_adj, g_from_hits, xi_active, g_undir);

    for (uint32_t h : hits)
    {
        if (!frimByteMapGet(g_from_hits, h) || !frimByteMapGet(xi_active, h))
            continue;
        g_undir[virtual_seed].push_back(h);
        g_undir[h].push_back(virtual_seed);
    }

    g_disc.clear();
    g_low.clear();
    g_parent.clear();
    g_subtree_has_root.clear();
    g_blocks_tarjan.clear();
    int timer = 0;
    frimPlatformTarjanDfsSparse(
        virtual_seed, g_undir, g_disc, g_low, g_parent,
        g_subtree_has_root, g_blocks_tarjan, sample.root, timer);

    for (uint32_t u : sample.nodes)
    {
        if (u == sample.root)
            continue;
        if (!frimByteMapGet(xi_active, u))
            continue;
        if (frimByteMapGet(g_from_hits, u) && frimByteMapGet(g_blocks_tarjan, u))
            blocks[u] = 1;
    }
    return true;
}

bool frimRRGraphSamplePlatformInRR(const FrimRRGraphSample& sample)
{
    return frimRRGraphSampleHasHits(sample);
}

void frimRRGraphFillBlocks(
    const FrimRRGraphSample& sample,
    FrimByteMap& blocks)
{
    blocks.clear();
    if (sample.nodes.empty() || !frimRRGraphSampleHasHits(sample))
        return;

    FrimByteMap node_active;
    for (uint32_t u : sample.nodes)
        node_active[u] = 1;

    frimRRGraphRunBlockPipeline(sample, blocks, node_active);
}

/// On a fixed xi≡1 RR-graph, flip one xi_lo coin per node, keep subgraph where
/// hit seeds and root stay connected, then Tarjan AP on the masked cone.
bool frimRRGraphFillBlocksWithXiLoMask(
    const FrimRRGraphSample& sample,
    double xi_lo,
    FrimByteMap& blocks)
{
    blocks.clear();
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return false;

    g_xi_active.clear();
    for (uint32_t u : sample.nodes)
    {
        if (u == sample.root)
            g_xi_active[u] = 1;
        else
            g_xi_active[u] = (dsfmt_gv_genrand_close_open() <= xi_lo) ? 1 : 0;
    }

    return frimRRGraphRunBlockPipeline(sample, blocks, g_xi_active);
}

void frimRRGraphFillXiActiveCrn(
    const FrimRRGraphSample& sample,
    const std::vector<double>& xi,
    FrimByteMap& xi_active,
    uint32_t override_u = UINT32_MAX,
    double override_xi = 0.0)
{
    xi_active.clear();
    for (uint32_t u : sample.nodes)
    {
        if (u == sample.root)
        {
            xi_active[u] = 1;
            continue;
        }
        const double xiu = (u == override_u)
            ? override_xi
            : ((u < xi.size()) ? xi[u] : 0.0);
        const auto gate_it = sample.xi_gate.find(u);
        if (gate_it != sample.xi_gate.end())
        {
            xi_active[u] =
                (static_cast<double>(gate_it->second) <= xiu) ? 1 : 0;
        }
        else
        {
            xi_active[u] = 1;
        }
    }
}

bool frimRRGraphSampleConnectedCrn(
    const FrimRRGraphSample& sample,
    const std::vector<double>& xi,
    uint32_t override_u = UINT32_MAX,
    double override_xi = 0.0)
{
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return false;

    frimRRGraphFillXiActiveCrn(sample, xi, g_xi_active, override_u, override_xi);
    g_from_hits.clear();
    frimFillReachFromMultiSource(
        sample.hit_nodes, sample.forward_adj, g_xi_active, g_from_hits);
    return frimByteMapGet(g_from_hits, sample.root);
}

/// Per-sample J contribution on fixed RR-graph: (1-lam_root*xi_root) if connected.
double frimRRGraphSampleJWeight(
    const FrimRRGraphSample& sample,
    const std::vector<double>& xi,
    const std::vector<double>& lam,
    uint32_t override_u = UINT32_MAX,
    double override_xi = 0.0)
{
    if (!frimRRGraphSampleConnectedCrn(sample, xi, override_u, override_xi))
        return 0.0;

    const uint32_t v = sample.root;
    if (v >= lam.size())
        return 0.0;

    const double xiv = (v == override_u)
        ? override_xi
        : ((v < xi.size()) ? xi[v] : 0.0);
    return 1.0 - lam[v] * xiv;
}


struct FrimApMonotoneCheck
{
    size_t hit_samples = 0;
    size_t masked_connected = 0;
    size_t violations = 0;
};

FrimApMonotoneCheck frimCheckApBlockMonotonicity(
    const std::vector<FrimRRGraphSample>& samples,
    double xi_lo)
{
    FrimApMonotoneCheck report;
    for (const auto& sample : samples)
    {
        if (!sample.hit || sample.nodes.empty())
            continue;
        report.hit_samples++;

        FrimByteMap full_blocks;
        frimRRGraphFillBlocks(sample, full_blocks);

        FrimByteMap masked_blocks;
        if (!frimRRGraphFillBlocksWithXiLoMask(sample, xi_lo, masked_blocks))
            continue;
        report.masked_connected++;

        for (uint32_t u : sample.nodes)
        {
            if (u == sample.root)
                continue;
            if (frimByteMapGet(full_blocks, u) && !frimByteMapGet(masked_blocks, u))
                report.violations++;
        }
    }
    return report;
}

/// xi≡xi_uniform on CRN gates: active iff xi_gate[u] <= xi_uniform (root always on).
static void frimRRGraphFillXiActiveUniformCrn(
    const FrimRRGraphSample& sample,
    double xi_uniform,
    FrimByteMap& xi_active)
{
    xi_active.clear();
    for (uint32_t u : sample.nodes)
    {
        if (u == sample.root)
        {
            xi_active[u] = 1;
            continue;
        }
        const auto gate_it = sample.xi_gate.find(u);
        if (gate_it != sample.xi_gate.end())
        {
            xi_active[u] =
                (static_cast<double>(gate_it->second) <= xi_uniform) ? 1 : 0;
        }
        else
        {
            xi_active[u] = 1;
        }
    }
}

/// CRN active mask: xi_v=xi_lo for v!=flip_u, xi_{flip_u}=1 (root always on).
/// Used to require q-hit when only one coordinate is flipped to 1 for K_lb.
static bool frimRRGraphSampleQHitXiLoExceptOne(
    const FrimRRGraphSample& sample,
    double xi_lo,
    uint32_t flip_u)
{
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return false;

    g_xi_active.clear();
    for (uint32_t u : sample.nodes)
    {
        if (u == sample.root)
        {
            g_xi_active[u] = 1;
            continue;
        }
        const double xiu = (u == flip_u) ? 1.0 : xi_lo;
        const auto gate_it = sample.xi_gate.find(u);
        if (gate_it != sample.xi_gate.end())
        {
            g_xi_active[u] =
                (static_cast<double>(gate_it->second) <= xiu) ? 1 : 0;
        }
        else
        {
            g_xi_active[u] = 1;
        }
    }

    g_from_hits.clear();
    frimFillReachFromMultiSource(
        sample.hit_nodes, sample.forward_adj, g_xi_active, g_from_hits);
    return frimByteMapGet(g_from_hits, sample.root);
}

static bool frimRRGraphSampleConnectedUniformCrn(
    const FrimRRGraphSample& sample,
    double xi_uniform)
{
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return false;

    frimRRGraphFillXiActiveUniformCrn(sample, xi_uniform, g_xi_active);
    g_from_hits.clear();
    frimFillReachFromMultiSource(
        sample.hit_nodes, sample.forward_adj, g_xi_active, g_from_hits);
    return frimByteMapGet(g_from_hits, sample.root);
}

struct FrimFixedXiLoSampleCtx
{
    bool connected = false;
};

/// Per sample: xi≡xi_lo connectivity + from_hits for pass-through (CRN xi_gate).
static void frimRRGraphPrepareFixedXiLo(
    const FrimRRGraphSample& sample,
    double xi_lo,
    FrimFixedXiLoSampleCtx& ctx)
{
    ctx.connected = false;
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return;

    frimRRGraphFillXiActiveUniformCrn(sample, xi_lo, g_xi_active_xilo);
    g_from_hits_xilo.clear();
    frimFillReachFromMultiSource(
        sample.hit_nodes, sample.forward_adj, g_xi_active_xilo, g_from_hits_xilo);
    ctx.connected = frimByteMapGet(g_from_hits_xilo, sample.root);
}

/// Tarjan blocks on xi≡1 (CRN xi_gate); used for K_lb only.
static void frimRRGraphPrepareBlocksXi1(const FrimRRGraphSample& sample)
{
    g_blocks_xi1.clear();
    if (!frimRRGraphSampleHasHits(sample) || sample.nodes.empty())
        return;

    frimRRGraphFillXiActiveUniformCrn(sample, 1.0, g_xi_active);
    frimRRGraphRunBlockPipeline(sample, g_blocks_xi1, g_xi_active);
}

static void frimAccumulateKLBPruneFromSamples(
    const std::vector<FrimRRGraphSample>& samples,
    const std::vector<double>& lam,
    double xi_lo,
    size_t num_v,
    std::vector<double>* hit_root_xi1,
    std::vector<double>* k_lb_accum,
    std::vector<double>* hit_root_xilo,
    std::vector<double>* k_ub_accum)
{
    FrimFixedXiLoSampleCtx ctx;
    for (const auto& sample : samples)
    {
        if (!sample.hit || sample.nodes.empty())
            continue;

        const uint32_t v = sample.root;
        if (v >= num_v || v >= lam.size())
            continue;

        if (hit_root_xi1
            && frimRRGraphSampleConnectedUniformCrn(sample, 1.0))
            (*hit_root_xi1)[v] += 1.0;

        frimRRGraphPrepareFixedXiLo(sample, xi_lo, ctx);
        if (hit_root_xilo && ctx.connected)
            (*hit_root_xilo)[v] += 1.0;

        // K as ∂J/∂xi_u: no xi_lo pass-through filter.
        // K_lb: Tarjan blocks at ξ≡1, count u only if q still hits when
        //       ξ_v=ξ_lo (∀v≠u) and ξ_u=1 (single-flip q-hit).
        // K_ub: all q-hit samples containing u (GRRG non-empty ⇒ q-hit at ξ=1),
        //       weight (1-λ_root·ξ_lo).
        if (k_lb_accum)
            frimRRGraphPrepareBlocksXi1(sample);

        const double coef_lb = 1.0 - lam[v];
        const double coef_ub = 1.0 - lam[v] * xi_lo;
        for (uint32_t u : sample.nodes)
        {
            if (u == v || u >= num_v)
                continue;

            if (k_lb_accum && frimByteMapGet(g_blocks_xi1, u)
                && frimRRGraphSampleQHitXiLoExceptOne(sample, xi_lo, u))
                (*k_lb_accum)[u] += coef_lb;
            if (k_ub_accum)
                (*k_ub_accum)[u] += coef_ub;
        }
    }
}

static void frimFinalizeKMinusLPrune(
    FrimPruneResult& out,
    const std::vector<double>& lam,
    const std::vector<double>& hit_root_xi1,
    double scale)
{
    out.K_lb = out.K;
    out.L_ub.assign(out.K.size(), 0.0);
    out.alpha_lb.assign(out.K.size(), 0);
    out.pruned_to_one.assign(out.K.size(), 0);
    out.num_pruned = 0;
    out.num_uncertain = 0;

    for (size_t u = 0; u < out.K.size(); u++)
    {
        out.L_ub[u] = lam[u] * hit_root_xi1[u] * scale;
        out.alpha_lb[u] = out.K_lb[u] - out.L_ub[u];
        if (out.alpha_lb[u] > 0.0)
        {
            out.pruned_to_one[u] = 1;
            out.num_pruned++;
        }
        else
        {
            out.num_uncertain++;
        }
    }
}

static void frimFinalizeKUbLlbPrune(
    FrimPruneLoResult& out,
    const std::vector<double>& lam,
    const std::vector<double>& hit_root_xilo,
    double scale)
{
    out.K_ub = out.K;
    out.L_lb.assign(out.K.size(), 0.0);
    out.alpha_ub.assign(out.K.size(), 0.0);
    out.pruned_to_lo.assign(out.K.size(), 0);
    out.num_pruned = 0;
    out.num_uncertain = 0;

    for (size_t u = 0; u < out.K.size(); u++)
    {
        out.L_lb[u] = lam[u] * hit_root_xilo[u] * scale;
        out.alpha_ub[u] = out.K_ub[u] - out.L_lb[u];
        if (out.alpha_ub[u] < 0.0)
        {
            out.pruned_to_lo[u] = 1;
            out.num_pruned++;
        }
        else
        {
            out.num_uncertain++;
        }
    }
}

static void frimLogConflictBoundSamples(
    const FrimPruneResult& hi,
    const FrimPruneLoResult& lo,
    size_t max_samples,
    const char* log_tag)
{
    size_t printed = 0;
    const size_t n = std::min(hi.K_lb.size(), lo.K_ub.size());
    for (size_t u = 0; u < n && printed < max_samples; u++)
    {
        const bool fix_one =
            (u < hi.pruned_to_one.size()) && hi.pruned_to_one[u];
        const bool fix_lo =
            (u < lo.pruned_to_lo.size()) && lo.pruned_to_lo[u];
        if (!fix_one || !fix_lo)
            continue;

        const double k_lb = hi.K_lb[u];
        const double l_ub = hi.L_ub[u];
        const double k_ub = lo.K_ub[u];
        const double l_lb = lo.L_lb[u];
        std::cout << log_tag << " conflict sample u=" << u
                  << " K_lb=" << k_lb
                  << " L_ub=" << l_ub
                  << " K_ub=" << k_ub
                  << " L_lb=" << l_lb
                  << " (alpha_lb=" << (k_lb - l_ub)
                  << ", alpha_ub=" << (k_ub - l_lb) << ")"
                  << std::endl;
        printed++;
    }
}

static void frimAnalyzeCombinedKLBFromSamples(
    const std::vector<FrimRRGraphSample>& samples,
    const std::vector<double>& lam,
    double xi_lo,
    size_t num_v,
    double scale,
    FrimPruneResult& hi,
    FrimPruneLoResult& lo)
{
    hi.num_rr = samples.size();
    lo.num_rr = samples.size();
    lo.xi_lo = xi_lo;
    hi.hoeffding_margin = 0.0;
    lo.hoeffding_margin = 0.0;

    std::vector<double> hit_root_xi1(num_v, 0.0);
    std::vector<double> hit_root_xilo(num_v, 0.0);
    hi.K.assign(num_v, 0.0);
    lo.K.assign(num_v, 0.0);

    frimAccumulateKLBPruneFromSamples(
        samples, lam, xi_lo, num_v,
        &hit_root_xi1, &hi.K, &hit_root_xilo, &lo.K);

    for (size_t u = 0; u < num_v; u++)
    {
        hi.K[u] *= scale;
        lo.K[u] *= scale;
    }

    frimFinalizeKMinusLPrune(hi, lam, hit_root_xi1, scale);
    frimFinalizeKUbLlbPrune(lo, lam, hit_root_xilo, scale);
}

} // namespace

bool Alg::frimRRGraphNodeBlocks(
    const FrimRRGraphSample& sample,
    uint32_t u)
{
    if (u == sample.root || sample.nodes.empty())
        return false;

    FrimByteMap blocks;
    frimRRGraphFillBlocks(sample, blocks);
    return frimByteMapGet(blocks, u);
}

static double frimNodeInfluenceProbUb(
    uint32_t u,
    const std::vector<double>& q,
    const Graph& reverse_graph)
{
    const double qu = (u < q.size()) ? q[u] : 0.0;
    double prod_p = 1.0;
    if (u < reverse_graph.size())
    {
        for (const auto& nbr : reverse_graph[u])
            prod_p *= static_cast<double>(nbr.second);
    }
    const double pr = qu + (1.0 - qu) * (1.0 - prod_p);
    return std::min(1.0, std::max(0.0, pr));
}

static double frimComputeHoeffdingMargin(
    double margin_coef,
    size_t num_samples,
    double hoeffding_delta,
    double margin_scale,
    double margin_override)
{
    if (margin_override >= 0.0)
        return margin_override;
    if (num_samples == 0 || hoeffding_delta <= 0.0)
        return 0.0;
    return margin_coef
        * std::sqrt(std::log(2.0 / hoeffding_delta)
                    / (2.0 * static_cast<double>(num_samples)))
        * margin_scale;
}

FrimPruneResult Alg::frimAnalyzeKMinusLPruneFromSamples(
    const std::vector<FrimRRGraphSample>& samples,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    double eps,
    double hoeffding_delta,
    double hoeffding_margin_scale,
    double hoeffding_margin_override) const
{
    FrimPruneResult out;
    out.num_rr = samples.size();

    if (_numV == 0 || samples.empty() || tau.size() != _numV
        || lam.size() != _numV || q.size() != _numV)
        return out;

    out.tau_sum = frimTauSum(tau);
    if (out.tau_sum <= 0.0)
        return out;

    const double scale = out.tau_sum / static_cast<double>(samples.size());
    (void)eps;
    (void)hoeffding_delta;
    (void)hoeffding_margin_scale;
    (void)hoeffding_margin_override;
    out.hoeffding_margin = 0.0;

    std::vector<double> hit_root_xi1(_numV, 0.0);
    out.K.assign(_numV, 0.0);

    frimAccumulateKLBPruneFromSamples(
        samples, lam, xi_lo, _numV,
        &hit_root_xi1, &out.K, nullptr, nullptr);

    for (size_t u = 0; u < _numV; u++)
        out.K[u] *= scale;

    frimFinalizeKMinusLPrune(out, lam, hit_root_xi1, scale);

    /*
    // --- previous bound (kept for reference) ---
    out.K.assign(_numV, 0.0);
    double max_hit_coef = 0.0;
    for (const auto& sample : samples)
    {
        if (!sample.hit || sample.nodes.empty())
            continue;
        const uint32_t v = sample.root;
        if (v >= lam.size())
            continue;
        const double coef = 1.0 - lam[v];
        max_hit_coef = std::max(max_hit_coef, coef);
        FrimByteMap blocks_old;
        frimRRGraphFillBlocks(sample, blocks_old);
        for (uint32_t u : sample.nodes)
        {
            if (u == v || u >= _numV)
                continue;
            if (frimByteMapGet(blocks_old, u))
                out.K[u] += coef;
        }
    }
    for (size_t u = 0; u < _numV; u++)
        out.K[u] *= scale;
    const double margin_coef = scale * max_hit_coef;
    (void)margin_coef;
    out.hoeffding_margin = eps;
    out.K_lb.assign(_numV, 0.0);
    out.L_ub.assign(_numV, 0.0);
    for (size_t u = 0; u < _numV; u++)
    {
        out.K_lb[u] = out.K[u] - out.hoeffding_margin;
        const double pr_ub =
            frimNodeInfluenceProbUb(static_cast<uint32_t>(u), q, _hyperGraph._graph);
        out.L_ub[u] = tau[u] * lam[u] * pr_ub;
        out.alpha_lb[u] = out.K_lb[u] - out.L_ub[u];
    }
    */

    (void)q;
    return out;
}

FrimPruneResult Alg::frimAnalyzeKMinusLPrune(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr,
    double eps,
    double hoeffding_delta) const
{
    if (_numV == 0 || num_rr == 0)
        return FrimPruneResult();

    const std::vector<FrimRRGraphSample> samples =
        frimBuildTauRootRRGraphSamplesAtXiOne(num_rr, tau, q);
    return frimAnalyzeKMinusLPruneFromSamples(
        samples, tau, lam, q, xi_lo, eps, hoeffding_delta);
}

void Alg::frimRunKMinusLPruneAnalysis(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr,
    double eps,
    double hoeffding_delta) const
{
    Timer timer;
    const FrimPruneResult result =
        frimAnalyzeKMinusLPrune(tau, lam, q, xi_lo, num_rr, eps, hoeffding_delta);

    size_t num_pointwise = 0;
    for (size_t u = 0; u < _numV; u++)
    {
        if (result.K[u] - result.L_ub[u] > 0.0)
            num_pointwise++;
    }

    std::cout << "[FRIM-PRUNE] xi ≡ 1, tau-root RR-graph, scale T/R"
              << ", R=" << num_rr << ", xi_lo=" << xi_lo << std::endl;
    std::cout << "[FRIM-PRUNE] L upper bound: lambda_u * (# q-hit RR-graphs with root u at xi=1)"
              << std::endl;
    std::cout << "[FRIM-PRUNE] K lower bound: Tarjan blocks at xi=1, count u only if"
              << " q still hits when xi_v=xi_lo (v!=u) and xi_u=1;"
              << " weight (1-lambda_root); no margin"
              << std::endl;
    std::cout << "[FRIM-PRUNE] prune if K_lb - L_ub > 0  =>  lock xi_u = 1"
              << std::endl;
    std::cout << "[FRIM-PRUNE] pruned: " << result.num_pruned
              << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE] uncertain (need search): " << result.num_uncertain
              << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE] pointwise K_hat - L_ub > 0: "
              << num_pointwise << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE] elapsed=" << timer.get_total_time() << "s" << std::endl;
}

FrimPruneLoResult Alg::frimAnalyzeKUbLlbAtXiLoFromSamples(
    const std::vector<FrimRRGraphSample>& samples,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    double eps,
    double hoeffding_delta,
    double hoeffding_margin_scale,
    double hoeffding_margin_override) const
{
    FrimPruneLoResult out;
    out.num_rr = samples.size();
    out.xi_lo = xi_lo;

    if (_numV == 0 || samples.empty() || tau.size() != _numV
        || lam.size() != _numV || q.size() != _numV)
        return out;

    out.tau_sum = frimTauSum(tau);
    if (out.tau_sum <= 0.0)
        return out;

    const double scale = out.tau_sum / static_cast<double>(samples.size());
    (void)eps;
    (void)hoeffding_delta;
    (void)hoeffding_margin_scale;
    (void)hoeffding_margin_override;
    out.hoeffding_margin = 0.0;

    std::vector<double> hit_root_xilo(_numV, 0.0);
    out.K.assign(_numV, 0.0);

    frimAccumulateKLBPruneFromSamples(
        samples, lam, xi_lo, _numV,
        nullptr, nullptr, &hit_root_xilo, &out.K);

    for (size_t u = 0; u < _numV; u++)
        out.K[u] *= scale;

    frimFinalizeKUbLlbPrune(out, lam, hit_root_xilo, scale);

    /*
    // --- previous bound (kept for reference) ---
    out.K.assign(_numV, 0.0);
    double max_hit_coef = 0.0;
    for (const auto& sample : samples)
    {
        if (!sample.hit || sample.nodes.empty())
            continue;
        const uint32_t v = sample.root;
        if (v >= lam.size())
            continue;
        FrimByteMap blocks_old;
        if (!frimRRGraphFillBlocksWithXiLoMask(sample, xi_lo, blocks_old))
            continue;
        const double coef = 1.0 - lam[v] * xi_lo;
        max_hit_coef = std::max(max_hit_coef, coef);
        for (uint32_t u : sample.nodes)
        {
            if (u == v || u >= _numV)
                continue;
            if (frimByteMapGet(blocks_old, u))
                out.K[u] += coef;
        }
    }
    for (size_t u = 0; u < _numV; u++)
        out.K[u] *= scale;
    out.hoeffding_margin = eps;
    out.K_ub.assign(_numV, 0.0);
    out.L_lb.assign(_numV, 0.0);
    for (size_t u = 0; u < _numV; u++)
    {
        out.K_ub[u] = out.K[u] + out.hoeffding_margin;
        const double qu = (u < q.size()) ? q[u] : 0.0;
        out.L_lb[u] = tau[u] * lam[u] * qu;
        out.alpha_ub[u] = out.K_ub[u] - out.L_lb[u];
    }
    */

    (void)q;
    return out;
}

FrimPruneLoResult Alg::frimAnalyzeKUbLlbAtXiLo(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr,
    double eps,
    double hoeffding_delta) const
{
    if (_numV == 0 || num_rr == 0)
        return FrimPruneLoResult();

    const std::vector<FrimRRGraphSample> samples =
        frimBuildTauRootRRGraphSamplesAtXiOne(num_rr, tau, q);
    return frimAnalyzeKUbLlbAtXiLoFromSamples(samples, tau, lam, q, xi_lo, eps, hoeffding_delta);
}

void Alg::frimRunKUbLlbPruneAnalysis(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr,
    double eps,
    double hoeffding_delta) const
{
    Timer timer;
    const FrimPruneLoResult result =
        frimAnalyzeKUbLlbAtXiLo(tau, lam, q, xi_lo, num_rr, eps, hoeffding_delta);

    size_t num_pointwise = 0;
    for (size_t u = 0; u < _numV; u++)
    {
        if (result.K[u] - result.L_lb[u] < 0.0)
            num_pointwise++;
    }

    std::cout << "[FRIM-PRUNE-LO] xi ≡ " << xi_lo
              << ", tau-root RR-graph at xi=1, scale T/R, R="
              << result.num_rr << std::endl;
    std::cout << "[FRIM-PRUNE-LO] L lower bound: lambda_u * (# q-hit RR-graphs with root u at xi=xi_lo)"
              << std::endl;
    std::cout << "[FRIM-PRUNE-LO] K upper bound: exclude pass-through-u samples at xi=xi_lo;"
              << " sum (1-lambda_root*xi_lo) over remaining samples with u in RR-graph; no margin"
              << std::endl;
    std::cout << "[FRIM-PRUNE-LO] prune if K_ub - L_lb < 0  =>  lock xi_u = " << xi_lo
              << std::endl;
    std::cout << "[FRIM-PRUNE-LO] pruned: " << result.num_pruned
              << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE-LO] uncertain (need search): " << result.num_uncertain
              << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE-LO] pointwise K_hat - L_lb < 0 (no Hoeffding): "
              << num_pointwise << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE-LO] elapsed=" << timer.get_total_time() << "s" << std::endl;
}

void Alg::frimRunCombinedPruneAnalysis(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    double xi_lo,
    size_t num_rr,
    int rand_seed,
    double eps,
    double hoeffding_delta) const
{
    Timer timer;
    const uint32_t seed = (rand_seed != 0)
        ? static_cast<uint32_t>(rand_seed)
        : static_cast<uint32_t>(time(nullptr));

    dsfmt_gv_init_gen_rand(seed);
    const std::vector<FrimRRGraphSample> samples =
        frimBuildTauRootRRGraphSamplesAtXiOne(num_rr, tau, q);
    const double tau_sum = frimTauSum(tau);
    const double scale = tau_sum / static_cast<double>(samples.size());

    FrimPruneResult hi;
    FrimPruneLoResult lo;
    hi.tau_sum = lo.tau_sum = tau_sum;
    frimAnalyzeCombinedKLBFromSamples(
        samples, lam, xi_lo, _numV, scale, hi, lo);

    size_t to_one = 0;
    size_t to_lo = 0;
    size_t conflict = 0;
    size_t uncertain = 0;

    for (size_t u = 0; u < _numV; u++)
    {
        const bool fix_one = (u < hi.pruned_to_one.size()) && hi.pruned_to_one[u];
        const bool fix_lo = (u < lo.pruned_to_lo.size()) && lo.pruned_to_lo[u];
        if (fix_one && fix_lo)
            conflict++;
        else if (fix_one)
            to_one++;
        else if (fix_lo)
            to_lo++;
        else
            uncertain++;
    }

    std::cout << "[FRIM-PRUNE-BOTH] R=" << num_rr << ", xi_lo=" << xi_lo << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] alone fix xi=1:  " << hi.num_pruned << " / " << _numV
              << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] alone fix xi=lo: " << lo.num_pruned << " / " << _numV
              << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] fix xi=1 (no conflict): " << to_one << " / " << _numV
              << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] fix xi=lo (no conflict): " << to_lo << " / " << _numV
              << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] conflict (both): " << conflict << " / " << _numV << std::endl;
    frimLogConflictBoundSamples(hi, lo, 5, "[FRIM-PRUNE-BOTH]");
    std::cout << "[FRIM-PRUNE-BOTH] uncertain: " << uncertain << " / " << _numV << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] total fixed (non-conflict): " << (to_one + to_lo)
              << " / " << _numV
              << " (" << (100.0 * static_cast<double>(to_one + to_lo) / static_cast<double>(_numV))
              << "%)" << std::endl;
    std::cout << "[FRIM-PRUNE-BOTH] elapsed=" << timer.get_total_time() << "s" << std::endl;
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

            const double xi_new = (result.xi[u] > 0.5) ? xi_lo : 1.0;
            if (std::abs(result.xi[u] - xi_new) <= 0.0)
                continue;

            std::vector<double> xi_flip = result.xi;
            xi_flip[u] = xi_new;
            const double J_stay = frimEstimateJ(result.xi, q, tau, lam, num_rr);
            const double J_flip = frimEstimateJ(xi_flip, q, tau, lam, num_rr);

            if (J_flip > J_stay + delta_tol)
            {
                result.xi[u] = xi_new;
                result.J_hat = J_flip;
                improved = true;
                updates_this_sweep++;
            }
        }

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
    result.J_method_rr_naive = result.J_hat;
    result.run_info.time_total_sec = timer.get_total_time();
    return result;
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

static const std::vector<uint32_t>& frimRRStructureHitList(
    const FrimRRStructureSample& sample)
{
    if (!sample.hit_nodes.empty())
        return sample.hit_nodes;
    static thread_local std::vector<uint32_t> legacy;
    legacy.clear();
    if (sample.hit)
        legacy.push_back(sample.hit_node);
    return legacy;
}

/// Forward propagation on the sampled RR subgraph: q-hit seeds -> root with CRN xi_gate.
static bool frimRRStructureSeedReachRoot(
    const FrimRRStructureSample& sample,
    const std::function<double(uint32_t)>& xiAt,
    FrimByteMap& reachable,
    std::unordered_map<uint32_t, uint32_t>* bfs_parent = nullptr)
{
    reachable.clear();
    if (bfs_parent)
        bfs_parent->clear();

    if (!sample.hit)
        return false;

    const std::vector<uint32_t>& hits = frimRRStructureHitList(sample);
    if (hits.empty())
        return false;

    std::queue<uint32_t> que;
    for (uint32_t h : hits)
    {
        if (!frimByteMapGet(reachable, h))
        {
            reachable[h] = 1;
            que.push(h);
        }
    }

    while (!que.empty())
    {
        const uint32_t u = que.front();
        que.pop();

        if (u != sample.root)
        {
            const auto gate_it = sample.xi_gate.find(u);
            if (gate_it == sample.xi_gate.end())
                continue;
            if (static_cast<double>(gate_it->second) > xiAt(u))
                continue;
        }

        const auto adj_it = sample.forward_adj.find(u);
        if (adj_it == sample.forward_adj.end())
            continue;

        for (uint32_t v : adj_it->second)
        {
            if (frimByteMapGet(reachable, v))
                continue;
            reachable[v] = 1;
            if (bfs_parent)
                (*bfs_parent)[v] = u;
            que.push(v);
        }
    }

    return frimByteMapGet(reachable, sample.root);
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

    FrimByteMap reachable;
    std::unordered_map<uint32_t, uint32_t> bfs_parent;
    if (!frimRRStructureSeedReachRoot(sample, xiAt, reachable, &bfs_parent))
        return 0.0;

    const std::vector<uint32_t>& hits = frimRRStructureHitList(sample);
    for (uint32_t hit_node : hits)
    {
        if (!frimByteMapGet(reachable, hit_node))
            continue;

        uint32_t cur = sample.root;
        if (!frimByteMapGet(reachable, cur))
            continue;

        double log_w = 0.0;
        while (cur != hit_node)
        {
            const auto pit = bfs_parent.find(cur);
            if (pit == bfs_parent.end())
            {
                log_w = 0.0;
                break;
            }
            const uint32_t p = pit->second;
            if (cur != sample.root && cur != hit_node)
            {
                const double xv = xiAt(cur);
                if (xv <= 0.0)
                    return 0.0;
                log_w += std::log(xv);
            }
            cur = p;
        }
        if (cur != hit_node)
            continue;

        if (sample.root >= tau.size() || sample.root >= lam.size())
            return 0.0;

        const double root_coef = 1.0 - lam[sample.root] * xiAt(sample.root);
        return root_coef * std::exp(log_w);
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

/// True when toggling xi[u] between xi_lo and 1 can change sample r's sweep weight.
/// Root always affects (1-lam_root*xi_root); non-root only when xi_gate > xi_lo.
static bool frimRRGraphSampleAffectsSweep(
    const FrimRRGraphSample& sample,
    uint32_t u,
    double xi_lo)
{
    if (u == sample.root)
        return true;
    const auto gate_it = sample.xi_gate.find(u);
    if (gate_it == sample.xi_gate.end())
        return false;
    return static_cast<double>(gate_it->second) > xi_lo + 1e-15;
}

/// frim_rr_graph: tau-root RR-graph + CRN xi_gate.
/// Per-sample weight = (1-lam_root*xi_root) if hit seeds connect to root, else 0.
/// J_hat = (tau_sum/R) * sum_r weight_r  (weighted coverage).
FrimXiResult Alg::frimMethodRRGraph(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol,
        bool combined_prune_init,
        double eps,
        double hoeffding_delta,
        double hoeffding_margin_scale,
        double hoeffding_margin_override,
        bool rr_graph_gate_sweep_index,
        bool store_hit_only,
        const std::vector<uint32_t>* subsim_fix_one)
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_rr == 0 || q.size() != _numV
        || tau.size() != _numV || lam.size() != _numV)
        return result;

    if (result.xi.size() != _numV)
        result.xi.assign(_numV, 1.0);

    const double tau_sum = frimTauSum(tau);
    if (tau_sum <= 0.0)
        return result;

    const char* log_tag =
        combined_prune_init ? "[FRIM-RR-GRAPH-PRUNE]" : "[FRIM-RR-GRAPH]";

    auto computeSampleWeight = [&](const FrimRRGraphSample& sample,
                                   const std::vector<double>& xi_local,
                                   uint32_t override_node = UINT32_MAX,
                                   double override_xi = 0.0) -> double
    {
        return frimRRGraphSampleJWeight(
            sample, xi_local, lam, override_node, override_xi);
    };

    Timer timer;
    std::vector<FrimRRGraphSample> samples;
    samples.reserve(num_rr);
    result.run_info.time_sample_sec = 0.0;

#if 0  // warm-start disabled
    size_t warm_n = 0;
    if (rr_graph_warm_start && num_rr > 1)
    {
        const size_t div = (rr_graph_warm_start_div > 0) ? rr_graph_warm_start_div : 10;
        warm_n = std::max<size_t>(1, num_rr / div);
        warm_n = std::min(warm_n, num_rr - 1);
    }

    if (warm_n > 0)
    {
        std::cout << log_tag << " warm-start: sampling first " << warm_n
                  << " tau-root RR-graphs..." << std::endl;
        Timer timer_warm_sample;
        samples = frimBuildTauRootRRGraphSamplesCrn(warm_n, tau, q);
        result.run_info.time_sample_sec += timer_warm_sample.get_operation_time();

        std::vector<std::vector<size_t>> warm_node_samples(_numV);
        for (size_t r = 0; r < samples.size(); r++)
        {
            for (uint32_t u : samples[r].nodes)
            {
                if (u < _numV)
                    warm_node_samples[u].push_back(r);
            }
        }

        std::vector<double> warm_weights(samples.size(), 0.0);
        double warm_cover_sum = 0.0;
        for (size_t r = 0; r < samples.size(); r++)
        {
            warm_weights[r] = computeSampleWeight(samples[r], result.xi);
            warm_cover_sum += warm_weights[r];
        }

        const double warm_scale = tau_sum / static_cast<double>(warm_n);
        result.J_hat = warm_scale * warm_cover_sum;

        auto warmRecomputeAffected = [&](uint32_t u, double xi_u_val) -> double
        {
            double affected = 0.0;
            for (size_t r : warm_node_samples[u])
                affected += computeSampleWeight(samples[r], result.xi, u, xi_u_val);
            return affected;
        };

        std::vector<uint32_t> warm_active(_numV);
        std::iota(warm_active.begin(), warm_active.end(), 0U);

        const size_t warm_sweep_count =
            (rr_graph_warm_start_sweeps > 0) ? rr_graph_warm_start_sweeps : 1;
        std::cout << log_tag << " warm-start: " << warm_sweep_count
                  << " sweep(s) on first " << warm_n
                  << " samples, J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        for (size_t warm_s = 0; warm_s < warm_sweep_count; warm_s++)
        {
            bool warm_improved = false;
            std::vector<uint32_t> warm_order = warm_active;
            for (size_t i = warm_order.size(); i > 1; i--)
            {
                const size_t j =
                    dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i));
                std::swap(warm_order[i - 1], warm_order[j]);
            }

            size_t warm_updates = 0;
            for (uint32_t u : warm_order)
            {
                double affected_old = 0.0;
                for (size_t r : warm_node_samples[u])
                    affected_old += warm_weights[r];

                const double affected_one = warmRecomputeAffected(u, 1.0);
                const double affected_lo = warmRecomputeAffected(u, xi_lo);

                const double J_one =
                    warm_scale * (warm_cover_sum - affected_old + affected_one);
                const double J_lo =
                    warm_scale * (warm_cover_sum - affected_old + affected_lo);

                const double xi_pick =
                    (J_one > J_lo + delta_tol) ? 1.0
                    : (J_lo > J_one + delta_tol) ? xi_lo
                    : (J_one >= J_lo ? 1.0 : xi_lo);

                if (std::abs(result.xi[u] - xi_pick) > 0.0)
                {
                    result.xi[u] = xi_pick;
                    const double affected_new =
                        (std::abs(xi_pick - 1.0) <= 0.0) ? affected_one : affected_lo;
                    warm_cover_sum = warm_cover_sum - affected_old + affected_new;
                    for (size_t r : warm_node_samples[u])
                        warm_weights[r] = computeSampleWeight(samples[r], result.xi);
                    result.J_hat = warm_scale * warm_cover_sum;
                    warm_improved = true;
                    warm_updates++;
                }
            }

            std::cout << log_tag << " warm-start sweep " << (warm_s + 1) << "/"
                      << warm_sweep_count << " done: updates=" << warm_updates
                      << ", improved=" << (warm_improved ? "yes" : "no")
                      << ", J_hat=" << result.J_hat
                      << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
        }

        const size_t rest_n = num_rr - warm_n;
        std::cout << log_tag << " warm-start: sampling remaining " << rest_n
                  << " tau-root RR-graphs..." << std::endl;
        Timer timer_rest_sample;
        const std::vector<FrimRRGraphSample> rest_samples =
            frimBuildTauRootRRGraphSamplesCrn(rest_n, tau, q);
        result.run_info.time_sample_sec += timer_rest_sample.get_operation_time();
        samples.insert(samples.end(), rest_samples.begin(), rest_samples.end());

        std::cout << log_tag << " sampling done (warm-start), total sample time="
                  << result.run_info.time_sample_sec << "s" << std::endl;
    }
    else
#endif
    {
        std::cout << log_tag << " sampling " << num_rr
                  << " tau-root RR-graphs (CRN xi_gate"
                  << (store_hit_only ? ", store_hit_only=on" : ", store_hit_only=off")
                  << ")..." << std::endl;
        size_t num_discarded = 0;
        samples = frimBuildTauRootRRGraphSamplesCrn(
            num_rr, tau, q, store_hit_only, &num_discarded);
        result.run_info.time_sample_sec = timer.get_operation_time();
        std::cout << log_tag << " sampling done, time="
                  << result.run_info.time_sample_sec << "s"
                  << ", stored=" << samples.size() << "/" << num_rr;
        if (store_hit_only && num_discarded > 0)
            std::cout << ", discarded_no_hit=" << num_discarded;
        std::cout << std::endl;
    }

    std::vector<uint8_t> xi_fixed(_numV, 0);
    std::vector<uint32_t> active_nodes;
    active_nodes.reserve(_numV);

    if (combined_prune_init)
    {
        Timer timer_prune;
        FrimPruneResult hi;
        FrimPruneLoResult lo;
        hi.tau_sum = lo.tau_sum = frimTauSum(tau);
        const double prune_scale =
            hi.tau_sum / static_cast<double>(num_rr);
        frimAnalyzeCombinedKLBFromSamples(
            samples, lam, xi_lo, _numV, prune_scale, hi, lo);

        size_t num_fix_one = 0;
        size_t num_fix_lo = 0;
        size_t num_conflict = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            const bool fix_one =
                (u < hi.pruned_to_one.size()) && hi.pruned_to_one[u];
            const bool fix_lo =
                (u < lo.pruned_to_lo.size()) && lo.pruned_to_lo[u];
            if (fix_one && fix_lo)
            {
                num_conflict++;
                continue;
            }
            if (fix_one)
            {
                result.xi[u] = 1.0;
                xi_fixed[u] = 1;
                num_fix_one++;
            }
            else if (fix_lo)
            {
                result.xi[u] = xi_lo;
                xi_fixed[u] = 1;
                num_fix_lo++;
            }
        }

        for (size_t u = 0; u < _numV; u++)
        {
            if (!xi_fixed[u])
                result.xi[u] = 1.0;
        }

        for (size_t u = 0; u < _numV; u++)
        {
            if (!xi_fixed[u])
                active_nodes.push_back(static_cast<uint32_t>(u));
        }

        std::cout << log_tag << " combined K/L prune (no margin)"
                  << ": fix xi=1=" << num_fix_one
                  << ", fix xi=lo=" << num_fix_lo
                  << ", conflict=" << num_conflict
                  << ", active=" << active_nodes.size() << "/" << _numV
                  << ", prune_time=" << timer_prune.get_operation_time() << "s"
                  << std::endl;
        frimLogConflictBoundSamples(hi, lo, 5, log_tag);
        size_t l_consistent = 0, l_inconsistent = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (u < hi.L_ub.size() && u < lo.L_lb.size())
            {
                if (hi.L_ub[u] + 1e-15 >= lo.L_lb[u])
                    l_consistent++;
                else
                    l_inconsistent++;
            }
        }
        std::cout << log_tag << " L_ub >= L_lb check: "
                  << l_consistent << " / " << (l_consistent + l_inconsistent)
                  << " (inconsistent=" << l_inconsistent << ")" << std::endl;
        size_t k_consistent = 0, k_inconsistent = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (u < hi.K_lb.size() && u < lo.K_ub.size())
            {
                if (lo.K_ub[u] + 1e-15 >= hi.K_lb[u])
                    k_consistent++;
                else
                    k_inconsistent++;
            }
        }
        std::cout << log_tag << " K_ub >= K_lb check: "
                  << k_consistent << " / " << (k_consistent + k_inconsistent)
                  << " (inconsistent=" << k_inconsistent << ")" << std::endl;
        size_t l_strict = 0, k_strict = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (u < hi.L_ub.size() && u < lo.L_lb.size()
                && hi.L_ub[u] > lo.L_lb[u] + 1e-15)
                l_strict++;
            if (u < hi.K_lb.size() && u < lo.K_ub.size()
                && lo.K_ub[u] > hi.K_lb[u] + 1e-15)
                k_strict++;
        }
        std::cout << log_tag << " L_ub > L_lb strict: " << l_strict << " / " << _numV << std::endl;
        std::cout << log_tag << " K_ub > K_lb strict: " << k_strict << " / " << _numV << std::endl;
        size_t alpha_lb_zero = 0, alpha_ub_zero = 0;
        size_t alpha_lb_pos = 0, alpha_lb_neg = 0;
        size_t alpha_ub_pos = 0, alpha_ub_neg = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (u >= hi.K_lb.size() || u >= hi.L_ub.size()
                || u >= lo.K_ub.size() || u >= lo.L_lb.size())
                continue;
            const double alb = hi.K_lb[u] - hi.L_ub[u];
            const double aub = lo.K_ub[u] - lo.L_lb[u];
            if (std::abs(alb) <= 1e-15)
                alpha_lb_zero++;
            else if (alb > 0.0)
                alpha_lb_pos++;
            else
                alpha_lb_neg++;
            if (std::abs(aub) <= 1e-15)
                alpha_ub_zero++;
            else if (aub > 0.0)
                alpha_ub_pos++;
            else
                alpha_ub_neg++;
        }
        std::cout << log_tag << " alpha_lb=K_lb-L_ub: >0 " << alpha_lb_pos
                  << ", =0 " << alpha_lb_zero
                  << ", <0 " << alpha_lb_neg << std::endl;
        std::cout << log_tag << " alpha_ub=K_ub-L_lb: >0 " << alpha_ub_pos
                  << ", =0 " << alpha_ub_zero
                  << ", <0 " << alpha_ub_neg << std::endl;

        // Conflict diagnostics: compare bound magnitudes on conflict nodes.
        size_t c_klb_gt_kub = 0, c_klb_eq_kub = 0, c_klb_lt_kub = 0;
        size_t c_lub_gt_llb = 0, c_lub_eq_llb = 0, c_lub_lt_llb = 0;
        size_t c_klb_gt_lub = 0, c_kub_lt_llb = 0;
        double sum_klb = 0, sum_kub = 0, sum_lub = 0, sum_llb = 0;
        size_t n_conflict_diag = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            const bool fix_one =
                (u < hi.pruned_to_one.size()) && hi.pruned_to_one[u];
            const bool fix_lo =
                (u < lo.pruned_to_lo.size()) && lo.pruned_to_lo[u];
            if (!fix_one || !fix_lo)
                continue;
            n_conflict_diag++;
            const double klb = hi.K_lb[u], kub = lo.K_ub[u];
            const double lub = hi.L_ub[u], llb = lo.L_lb[u];
            sum_klb += klb; sum_kub += kub; sum_lub += lub; sum_llb += llb;
            if (klb > kub + 1e-15) c_klb_gt_kub++;
            else if (std::abs(klb - kub) <= 1e-15) c_klb_eq_kub++;
            else c_klb_lt_kub++;
            if (lub > llb + 1e-15) c_lub_gt_llb++;
            else if (std::abs(lub - llb) <= 1e-15) c_lub_eq_llb++;
            else c_lub_lt_llb++;
            if (klb > lub) c_klb_gt_lub++;
            if (kub < llb) c_kub_lt_llb++;
        }
        if (n_conflict_diag > 0)
        {
            std::cout << log_tag << " conflict bound sizes (n=" << n_conflict_diag << "):"
                      << " mean K_lb=" << (sum_klb / n_conflict_diag)
                      << " K_ub=" << (sum_kub / n_conflict_diag)
                      << " L_ub=" << (sum_lub / n_conflict_diag)
                      << " L_lb=" << (sum_llb / n_conflict_diag)
                      << std::endl;
            std::cout << log_tag << " conflict K_lb ? K_ub: > "
                      << c_klb_gt_kub << ", = " << c_klb_eq_kub
                      << ", < " << c_klb_lt_kub << std::endl;
            std::cout << log_tag << " conflict L_ub ? L_lb: > "
                      << c_lub_gt_llb << ", = " << c_lub_eq_llb
                      << ", < " << c_lub_lt_llb << std::endl;
        }

        // Global: how often K_lb > K_ub (expected if different eval points)
        size_t g_klb_gt_kub = 0, g_klb_eq_kub = 0, g_klb_lt_kub = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (u >= hi.K_lb.size() || u >= lo.K_ub.size())
                continue;
            const double klb = hi.K_lb[u], kub = lo.K_ub[u];
            if (klb > kub + 1e-15) g_klb_gt_kub++;
            else if (std::abs(klb - kub) <= 1e-15) g_klb_eq_kub++;
            else g_klb_lt_kub++;
        }
        std::cout << log_tag << " global K_lb ? K_ub: > "
                  << g_klb_gt_kub << ", = " << g_klb_eq_kub
                  << ", < " << g_klb_lt_kub << std::endl;
    }
    else if (subsim_fix_one != nullptr && !subsim_fix_one->empty())
    {
        std::unordered_set<uint32_t> seed_set(
            subsim_fix_one->begin(), subsim_fix_one->end());
        size_t num_fix_one = 0;
        size_t num_fix_lo = 0;
        for (size_t u = 0; u < _numV; u++)
        {
            if (seed_set.count(static_cast<uint32_t>(u)))
            {
                result.xi[u] = 1.0;
                xi_fixed[u] = 1;
                num_fix_one++;
            }
            else
            {
                result.xi[u] = xi_lo;
                xi_fixed[u] = 1;
                num_fix_lo++;
            }
        }
        std::cout << log_tag << " subsim init: fix xi=1 on " << num_fix_one
                  << " seeds, fix xi=" << xi_lo << " on " << num_fix_lo
                  << " other nodes (no sweep)" << std::endl;
    }
    else
    {
        active_nodes.resize(_numV);
        std::iota(active_nodes.begin(), active_nodes.end(), 0U);
    }

    std::vector<std::vector<size_t>> node_samples(_numV);
    size_t index_pairs_full = 0;
    size_t index_pairs_kept = 0;
    for (size_t r = 0; r < samples.size(); r++)
    {
        const FrimRRGraphSample& sample = samples[r];
        for (uint32_t u : sample.nodes)
        {
            if (u >= _numV)
                continue;
            index_pairs_full++;
            if (rr_graph_gate_sweep_index
                && !frimRRGraphSampleAffectsSweep(sample, u, xi_lo))
                continue;
            node_samples[u].push_back(r);
            index_pairs_kept++;
        }
    }
    if (rr_graph_gate_sweep_index && index_pairs_full > 0)
    {
        const size_t skipped = index_pairs_full - index_pairs_kept;
        std::cout << log_tag << " gate-sweep-index: kept " << index_pairs_kept
                  << "/" << index_pairs_full << " (u,r) pairs, skipped "
                  << skipped << " (" << (100.0 * skipped / index_pairs_full)
                  << "%)" << std::endl;
    }

    auto computeSampleWeightMain = computeSampleWeight;

    // Weighted coverage: sum of per-RR-set weights (stored samples only)
    std::vector<double> sample_weights(samples.size(), 0.0);
    double weighted_cover_sum = 0.0;
    for (size_t r = 0; r < samples.size(); r++)
    {
        sample_weights[r] = computeSampleWeightMain(samples[r], result.xi);
        weighted_cover_sum += sample_weights[r];
    }

    const double scale = tau_sum / static_cast<double>(num_rr);
    result.J_hat = scale * weighted_cover_sum;
    std::cout << log_tag << " initial J_hat=" << result.J_hat
              << " (tau_sum/R * weighted_cover, R=" << num_rr
              << ", weighted_cover=" << weighted_cover_sum << ")"
              << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

    std::vector<uint32_t> sweep_nodes;
    sweep_nodes.reserve(active_nodes.size());
    for (uint32_t u : active_nodes)
    {
        if (!xi_fixed[u])
            sweep_nodes.push_back(u);
    }

    auto recomputeAffected = [&](uint32_t u, double xi_u_val) -> double
    {
        double affected = 0.0;
        for (size_t r : node_samples[u])
            affected += computeSampleWeightMain(samples[r], result.xi, u, xi_u_val);
        return affected;
    };

    size_t sample_update_count = 0;
    double time_sample_update_sec = 0.0;

    int sweeps_completed = 0;
    for (int sweep = 0; sweep < max_sweeps; sweep++)
    {
        sweeps_completed++;
        bool improved = false;
        std::vector<uint32_t> order = sweep_nodes;
        for (size_t i = order.size(); i > 1; i--)
        {
            const size_t j =
                dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i));
            std::swap(order[i - 1], order[j]);
        }

        std::cout << log_tag << " sweep " << (sweep + 1) << "/" << max_sweeps
                  << " start, active_nodes=" << order.size()
                  << ", J_hat=" << result.J_hat
                  << ", elapsed=" << timer.get_total_time() << "s" << std::endl;

        uint32_t currProgress = 0;
        size_t node_idx = 0;
        size_t updates_this_sweep = 0;
        const size_t active_n = order.size();
        for (uint32_t u : order)
        {
            node_idx++;
            if (active_n > 0 && node_idx * 100 >= active_n * currProgress)
            {
                std::cout << log_tag << " sweep " << (sweep + 1) << "/" << max_sweeps
                          << " nodes " << currProgress << "%, updates="
                          << updates_this_sweep
                          << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
                currProgress += 20;
            }

            if (xi_fixed[u])
                continue;

            double affected_old = 0.0;
            for (size_t r : node_samples[u])
                affected_old += sample_weights[r];

            const double affected_one = recomputeAffected(u, 1.0);
            const double affected_lo = recomputeAffected(u, xi_lo);

            const double J_one = scale * (weighted_cover_sum - affected_old + affected_one);
            const double J_lo = scale * (weighted_cover_sum - affected_old + affected_lo);

            const double xi_pick =
                (J_one > J_lo + delta_tol) ? 1.0
                : (J_lo > J_one + delta_tol) ? xi_lo
                : (J_one >= J_lo ? 1.0 : xi_lo);

            if (std::abs(result.xi[u] - xi_pick) > 0.0)
            {
                result.xi[u] = xi_pick;
                const double affected_new =
                    (std::abs(xi_pick - 1.0) <= 0.0) ? affected_one : affected_lo;
                weighted_cover_sum = weighted_cover_sum - affected_old + affected_new;
                Timer timer_sample_update;
                for (size_t r : node_samples[u])
                {
                    sample_weights[r] = computeSampleWeightMain(samples[r], result.xi);
                    sample_update_count++;
                }
                time_sample_update_sec += timer_sample_update.get_operation_time();
                result.J_hat = scale * weighted_cover_sum;
                improved = true;
                updates_this_sweep++;
            }
        }

        std::cout << log_tag << " sweep " << (sweep + 1) << "/" << max_sweeps
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
    result.J_method_rr_graph = scale * weighted_cover_sum;
    result.J_hat = result.J_method_rr_graph;
    result.run_info.time_total_sec = timer.get_total_time();
    if (sample_update_count > 0)
    {
        const double ms_per_100 =
            time_sample_update_sec / static_cast<double>(sample_update_count) * 100.0
            * 1000.0;
        std::cout << log_tag << " sample_update_stats: count="
                  << sample_update_count
                  << ", time=" << time_sample_update_sec << "s"
                  << ", ms_per_100_updates=" << ms_per_100 << std::endl;
    }
    return result;
}

FrimXiResult Alg::frim_solve_rr_graph(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol,
    bool rr_graph_gate_sweep_index,
    bool store_hit_only)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-RR-GRAPH] " << num_rr
              << " tau-root RR-graph samples, CRN {xi_lo,1} via xi_gate..."
              << (rr_graph_gate_sweep_index ? ", gate-sweep-index=on" : ", gate-sweep-index=off")
              << (store_hit_only ? ", store_hit_only=on" : ", store_hit_only=off")
              << std::endl;
    const FrimXiResult method_graph = frimMethodRRGraph(
        tau, lam, q, init, xi_lo, num_rr, max_sweeps, delta_tol,
        false, 0.1, 0.05, 1.0, -1.0,
        rr_graph_gate_sweep_index, store_hit_only, nullptr);
    result = method_graph;
    result.J_hat = result.J_method_rr_graph;
    std::cout << "[FRIM-RR-GRAPH] done: J_hat=" << result.J_method_rr_graph << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-RR-GRAPH] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

FrimXiResult Alg::frim_solve_subsim(
    double xi_lo,
    double subsim_eps,
    double subsim_delta)
{
    FrimXiResult result;
    if (_numV == 0)
        return result;

    const int seed_k = std::max(1, static_cast<int>(_numV / 100));
    const double delta = (subsim_delta > 0.0)
        ? subsim_delta
        : (1.0 / static_cast<double>(_numV));

    Timer timer("FRIM-SUBSIM");
    std::cout << "[FRIM-SUBSIM] SubSim IM seeds: k=n/100=" << seed_k
              << ", eps=" << subsim_eps << ", delta=" << delta << std::endl;
    subsimOnly(seed_k, subsim_eps, delta);
    RefreshHypergraph();

    std::cout << "[FRIM-SUBSIM] selected " << _vecSeed.size() << " seeds:" << std::endl;
    for (uint32_t u : _vecSeed)
        std::cout << "  seed node " << u << std::endl;

    result.xi.assign(_numV, xi_lo);
    for (uint32_t u : _vecSeed)
    {
        if (u < _numV)
            result.xi[u] = 1.0;
    }

    result.run_info.time_total_sec = timer.get_total_time();
    result.run_info.time_sample_sec = result.run_info.time_total_sec;
    result.run_info.time_solve_sec = 0.0;
    result.run_info.time_estimate_sec = 0.0;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-SUBSIM] heuristic xi: " << result.run_info.num_xi_one
              << " at 1, " << result.run_info.num_xi_lo << " at " << xi_lo
              << " (no FRIM estimate)" << std::endl;
    return result;
}

FrimXiResult Alg::frim_solve_rr_graph_prune(
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    const std::vector<double>& q,
    const std::vector<double>& xi_init,
    double xi_lo,
    size_t num_rr,
    int max_sweeps,
    double delta_tol,
    double eps,
    double hoeffding_delta,
    double hoeffding_margin_scale,
    double hoeffding_margin_override,
    bool rr_graph_gate_sweep_index,
    bool store_hit_only)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-RR-GRAPH-PRUNE] combined K/L prune init + J sweep, R="
              << num_rr << ", max_sweeps=" << max_sweeps
              << ", hoeffding_delta=" << hoeffding_delta;
    if (hoeffding_margin_override >= 0.0)
        std::cout << ", hoeffding_margin=" << hoeffding_margin_override << " (fixed)";
    else if (std::abs(hoeffding_margin_scale - 1.0) >= 1e-12)
        std::cout << ", hoeffding_margin_scale=" << hoeffding_margin_scale;
    std::cout << ", gate-sweep-index="
              << (rr_graph_gate_sweep_index ? "on" : "off")
              << (store_hit_only ? ", store_hit_only=on" : ", store_hit_only=off")
              << std::endl;

    result = frimMethodRRGraph(
        tau, lam, q, init, xi_lo, num_rr, max_sweeps, delta_tol,
        true, eps, hoeffding_delta, hoeffding_margin_scale, hoeffding_margin_override,
        rr_graph_gate_sweep_index, store_hit_only, nullptr);

    std::cout << "[FRIM-RR-GRAPH-PRUNE] done: J_hat=" << result.J_method_rr_graph << std::endl;
    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-RR-GRAPH-PRUNE] xi: " << result.run_info.num_xi_one
              << " at 1, " << result.run_info.num_xi_lo << " at " << xi_lo << std::endl;
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
    double delta_tol,
    CascadeModel model) const
{
    FrimXiResult result;
    result.xi = xi_init;

    if (_numV == 0 || num_mc == 0 || q.size() != _numV
        || tau.size() != _numV || lam.size() != _numV)
        return result;

    Timer timer;
    result.J_hat = frimEstimateBenefitJ(
        forwardGraph, result.xi, q, tau, lam, model, num_mc);
    result.run_info.time_sample_sec = timer.get_operation_time();
    std::cout << "[FRIM-MC-naive] initial J_hat=" << result.J_hat
              << " (benefit_inf_eval, mc=" << num_mc << ")"
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

        std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
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
                std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
                          << " nodes " << currProgress << "%, updates=" << updates_this_sweep
                          << ", elapsed=" << timer.get_total_time() << "s" << std::endl;
                currProgress += 20;
            }

            const double xi_new = (result.xi[u] > 0.5) ? xi_lo : 1.0;
            if (std::abs(result.xi[u] - xi_new) <= 0.0)
                continue;

            std::vector<double> xi_flip = result.xi;
            xi_flip[u] = xi_new;
            const double J_stay = frimEstimateBenefitJ(
                forwardGraph, result.xi, q, tau, lam, model, num_mc);
            const double J_flip = frimEstimateBenefitJ(
                forwardGraph, xi_flip, q, tau, lam, model, num_mc);

            if (J_flip > J_stay + delta_tol)
            {
                result.xi[u] = xi_new;
                result.J_hat = J_flip;
                improved = true;
                updates_this_sweep++;
            }
        }

        std::cout << "[FRIM-MC-naive] sweep " << (sweep + 1) << "/" << max_sweeps
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
    result.J_method_mc_naive = result.J_hat;
    result.run_info.time_total_sec = timer.get_total_time();
    return result;
}

double Alg::frimEstimateBenefitJ(
    const Graph& forwardGraph,
    const std::vector<double>& xi,
    const std::vector<double>& q,
    const std::vector<double>& tau,
    const std::vector<double>& lam,
    CascadeModel model,
    size_t num_mc) const
{
    return GraphBase::benefit_inf_eval(
        forwardGraph, q, tau, lam, model, xi,
        static_cast<uint32_t>(num_mc), false);
}

FrimXiResult Alg::frim_solve_mc_crn(
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

    std::cout << "[FRIM-MC-CRN] coordinate ascent with " << num_mc
              << " CRN live-edge samples (random seed)..." << std::endl;
    const FrimXiResult method_mc = frimMethodMC(
        forwardGraph, tau, lam, q, init, xi_lo, num_mc, max_sweeps, delta_tol);
    result = method_mc;
    result.J_hat = result.J_method_mc;
    std::cout << "[FRIM-MC-CRN] done: J_hat=" << result.J_method_mc << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-MC-CRN] xi: " << result.run_info.num_xi_one << " nodes at 1, "
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
    double delta_tol,
    CascadeModel model)
{
    FrimXiResult result;
    if (_numV == 0 || tau.size() != _numV || lam.size() != _numV || q.size() != _numV)
        return result;

    std::vector<double> init = xi_init;
    if (init.size() != _numV)
        init.assign(_numV, 1.0);

    std::cout << "[FRIM-MC-naive] benefit_inf_eval, fresh MC per J call (mc=" << num_mc
              << ", model=" << (model == LT ? "LT" : "IC") << ")..." << std::endl;
    const FrimXiResult method_naive = frimMethodMCNaive(
        forwardGraph, tau, lam, q, init, xi_lo, num_mc, max_sweeps, delta_tol, model);
    result = method_naive;
    result.J_hat = result.J_method_mc_naive;
    std::cout << "[FRIM-MC-naive] done: J_hat=" << result.J_method_mc_naive << std::endl;

    frimFillXiStats(result, xi_lo);
    std::cout << "[FRIM-MC-naive] xi: " << result.run_info.num_xi_one << " nodes at 1, "
              << result.run_info.num_xi_lo << " nodes at " << xi_lo << std::endl;
    return result;
}

// ===== end FRIM xi selection =====
