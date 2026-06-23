#pragma once

class HyperGraph
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    std::vector<bool> _vecVisitBool;
    Nodelist _vecVisitNode;
    std::vector<uint32_t> _vecOutdeg;
    std::vector<size_t> _cumDeg;
    std::vector<size_t> _cumInDeg;
    std::vector<double> _vecWeightedIndeg;
    std::vector<double> _cumWeightedInDeg;
    double _totalWeightedIndeg = 0.0;
    std::vector<double> _vecWeightedOutdeg;
    std::vector<double> _cumWeightedOutdeg;
    double _totalWeightedOutdeg = 0.0;
    double _wrrNormConst = 0.0;
    
    size_t _hit = 0;
    size_t _frimNumRR = 0;
    std::vector<double> _frimRRweights;
    double _numSamplesEval = 0;
    double _hyperedgeAvgEval = 0.0;

    bool _isVanilla = false;
    WRRSampleMode _wrrSampleMode = WRR_INDEG;
    bool _wrrBuildActive = false;

    void updateWrrNormConst()
    {
        switch (_wrrSampleMode) {
        case WRR_WINDEG:
            _wrrNormConst = _totalWeightedIndeg > 0.0 ? _totalWeightedIndeg : static_cast<double>(_numV);
            break;
        case WRR_WOUTDEG:
            _wrrNormConst = _totalWeightedOutdeg + static_cast<double>(_numV);
            break;
        case WRR_OUTDEG:
            _wrrNormConst = static_cast<double>(_numV + _numE);
            break;
        case WRR_UNIFORM:
            _wrrNormConst = static_cast<double>(_numV);
            break;
        case WRR_INDEG:
        default:
            _wrrNormConst = _numE > 0 ? static_cast<double>(_numE) : static_cast<double>(_numV);
            break;
        }
    }

    /// Initialization
    void InitHypergraph()
    {
        _numV = (uint32_t)_graph.size();
        _vecOutdeg = std::vector<uint32_t>(_numV, 0);
        _nodeGain = std::vector<double>(_numV + 1, 0.0);
        _numE = 0;
        for (auto& nbrs : _graph) {
            _numE += nbrs.size();
            for (auto& nbr : nbrs) {
                _vecOutdeg[nbr.first]++;
            }
        }
        _cumDeg = std::vector<size_t>(_numV, 0); _cumDeg[0] = _vecOutdeg[0];
        _cumInDeg = std::vector<size_t>(_numV, 0); _cumInDeg[0] = _graph[0].size();
        for (int i = 1; i < _numV; i++) {
            _cumDeg[i] = _cumDeg[i-1] + _vecOutdeg[i];
            _cumInDeg[i] = _cumInDeg[i-1] + _graph[i].size();
        }
        _vecWeightedIndeg = std::vector<double>(_numV, 0.0);
        _cumWeightedInDeg = std::vector<double>(_numV, 0.0);
        _totalWeightedIndeg = 0.0;
        for (int i = 0; i < _numV; i++) {
            for (auto& nbr : _graph[i]) {
                _vecWeightedIndeg[i] += nbr.second;
            }
            _totalWeightedIndeg += _vecWeightedIndeg[i];
        }
        if (_numV > 0) {
            _cumWeightedInDeg[0] = _vecWeightedIndeg[0];
            for (int i = 1; i < _numV; i++) {
                _cumWeightedInDeg[i] = _cumWeightedInDeg[i - 1] + _vecWeightedIndeg[i];
            }
        }
        _vecWeightedOutdeg = std::vector<double>(_numV, 0.0);
        _cumWeightedOutdeg = std::vector<double>(_numV, 0.0);
        _totalWeightedOutdeg = 0.0;
        for (int v = 0; v < _numV; v++) {
            for (auto& nbr : _graph[v]) {
                _vecWeightedOutdeg[nbr.first] += nbr.second;
            }
        }
        for (int i = 0; i < _numV; i++) {
            _totalWeightedOutdeg += _vecWeightedOutdeg[i];
        }
        if (_numV > 0) {
            _cumWeightedOutdeg[0] = _vecWeightedOutdeg[0] + 1.0;
            for (int i = 1; i < _numV; i++) {
                _cumWeightedOutdeg[i] = _cumWeightedOutdeg[i - 1] + _vecWeightedOutdeg[i] + 1.0;
            }
        }
        // revise the sampling probability
        for (int i = 0; i < _numV; i++) {
            _cumDeg[i] += (i+1);
        }
        assert(_cumDeg[_numV-1] == (_numV + _numE));
        assert(_cumInDeg[_numV-1] == _numE);
        updateWrrNormConst();
        _FRsets = FRsets(_numV + 1);
        _vecVisitBool = std::vector<bool>(_numV + 1);
        _vecVisitNode = Nodelist(_numV);
    }

public:
    /// _graph: reverse graph
    const Graph& _graph;
    /// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
    FRsets _FRsets;
    /// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
    RRsets _RRsets;
    std::vector<double> _RRweights;
    std::vector<double> _nodeGain;
    ProbDist _probDist = WEIGHTS;

    explicit HyperGraph(const Graph& graph) : _graph(graph)
    {
        InitHypergraph();
    }

    /// Set cascade model
    void set_prob_dist(const ProbDist dist)
    {
        _probDist = dist;
    }

    void set_vanilla_sample(const bool isVanilla)
    {
        _isVanilla = isVanilla;
    }

    void set_wrr_sample_mode(const WRRSampleMode mode)
    {
        _wrrSampleMode = mode;
        updateWrrNormConst();
    }

    WRRSampleMode get_wrr_sample_mode() const
    {
        return _wrrSampleMode;
    }

    double get_wrr_norm_const() const
    {
        return _wrrNormConst;
    }

    /// Returns the number of nodes in the graph.
    uint32_t get_nodes() const
    {
        return _numV;
    }

    /// Returns the number of edges in the graph.
    size_t get_edges() const
    {
        return _numE;
    }

    /// Returns the number of RR sets in the graph.
    size_t get_RR_sets_size() const
    {
        return _numRRsets;
    }

    uint32_t sampleNodeByOutDegree() {
        // Handle edge cases
        if (_numV == 0) return -1;
        
        // Check if all out-degrees are 0
        // bool allZero = std::all_of(_vecOutdeg.begin(), _vecOutdeg.end(), 
        //                         [](int deg) { return deg == 0; });
        // if (allZero) return -1;
        
        // 1. Compute cumulative distribution
        // done in init
        
        // 2. Generate random number in [0, totalEdges)
        double U = dsfmt_gv_genrand_close_open();  // [0, 1)
        double X = U * (_numV + _numE);                    // [0, totalEdges)
        
        if (X < _cumDeg[0]) return 0;
        // 3. Binary search to find the node
        // Find first index where cum[i] > X
        
        int left = 0, right = _numV-1;
        while (left < right) {
            int mid = left + (right - left) / 2;
            if (_cumDeg[mid] <= X) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        return right;  // Node index
    }

    uint32_t sampleNodeByInDegree() {
        // Handle edge cases
        if (_numV == 0) return -1;
        
        // Check if all in-degrees are 0
        // bool allZero = std::all_of(_vecIndeg.begin(), _vecIndeg.end(), 
        //                         [](int deg) { return deg == 0; });
        // if (allZero) return -1;
        
        // 1. Compute cumulative distribution
        // done in init
        
        // 2. Generate random number in [0, totalEdges)
        double U = dsfmt_gv_genrand_close_open();  // [0, 1)
        double X = U * _numE;                    // [0, totalEdges)
        
        // 3. Binary search to find the node
        // Find first index where cum[i] > X
        if (X < _cumInDeg[0]) return 0;
        int left = 0, right = _numV-1;
        while (left < right) {
            int mid = left + (right - left) / 2;
            if (_cumInDeg[mid] <= X) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        return right;  // Node index
    }

    uint32_t sampleNodeByWeightedInDegree() {
        if (_numV == 0) return -1;
        if (_totalWeightedIndeg <= 0.0) {
            return dsfmt_gv_genrand_uint32_range(_numV);
        }

        const double U = dsfmt_gv_genrand_close_open();
        const double X = U * _totalWeightedIndeg;

        if (X < _cumWeightedInDeg[0]) return 0;
        int left = 0, right = static_cast<int>(_numV) - 1;
        while (left < right) {
            const int mid = left + (right - left) / 2;
            if (_cumWeightedInDeg[mid] <= X) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        return static_cast<uint32_t>(right);
    }

    uint32_t sampleNodeByWeightedOutDegree() {
        if (_numV == 0) return -1;

        const double totalMass = _totalWeightedOutdeg + static_cast<double>(_numV);
        const double U = dsfmt_gv_genrand_close_open();
        const double X = U * totalMass;

        if (X < _cumWeightedOutdeg[0]) return 0;
        int left = 0, right = static_cast<int>(_numV) - 1;
        while (left < right) {
            const int mid = left + (right - left) / 2;
            if (_cumWeightedOutdeg[mid] <= X) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        return static_cast<uint32_t>(right);
    }

    uint32_t sampleWRRRoot() {
        switch (_wrrSampleMode) {
        case WRR_WINDEG:
            return sampleNodeByWeightedInDegree();
        case WRR_WOUTDEG:
            return sampleNodeByWeightedOutDegree();
        case WRR_OUTDEG:
            return sampleNodeByOutDegree();
        case WRR_UNIFORM:
            return dsfmt_gv_genrand_uint32_range(_numV);
        case WRR_INDEG:
        default:
            if (_numE == 0) {
                return dsfmt_gv_genrand_uint32_range(_numV);
            }
            return sampleNodeByInDegree();
        }
    }

    double wrrRootWeight(uint32_t uStart) const {
        switch (_wrrSampleMode) {
        case WRR_WINDEG:
            assert(_vecWeightedIndeg[uStart] > 0.0);
            return 1.0 / _vecWeightedIndeg[uStart];
        case WRR_WOUTDEG:
            return 1.0 / (_vecWeightedOutdeg[uStart] + 1.0);
        case WRR_OUTDEG:
            return 1.0 / (_vecOutdeg[uStart] + 1);
        case WRR_UNIFORM:
            return 1.0;
        case WRR_INDEG:
        default:
            assert(_graph[uStart].size() > 0);
            return 1.0 / _graph[uStart].size();
        }
    }

    template<typename BuildFn>
    void buildWRRBatch(const size_t numSamples, const size_t prevSize, BuildFn buildFn)
    {
        for (auto i = prevSize; i < numSamples; i++) {
            const uint32_t u = sampleWRRRoot();
            buildFn(u, i);
            this->_RRweights.push_back(wrrRootWeight(u));
        }
    }

    /// Generate a set of n RR sets
    void BuildRRsets(const size_t numSamples)
    {
        void (*func)(const uint32_t uStart, const size_t hyperIdx);

        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR set by vanilla method" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }

        if (_probDist == WC)
        {
            // std::cout << "Sample RR sets in WC model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetWeighted(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == UNIFORM)
        {
            // std::cout << "Sample RR sets in uniform model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == SKEWED || _probDist == WEIGHTS)
        {
            // std::cout << "Sample RR sets in skewed or weights case" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetSkewed(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else
        {
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);;
            }
        }

        return ;
    }

    void BuildWRRsets(const size_t numSamples)
    {
        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;
        _wrrBuildActive = true;

        if (_isVanilla)
        {
            buildWRRBatch(numSamples, prevSize, [this](uint32_t u, size_t i) {
                BuildOneRRset(u, i);
            });
        }
        else if (_probDist == WC)
        {
            buildWRRBatch(numSamples, prevSize, [this](uint32_t u, size_t i) {
                BuildOneRRsetWeighted(u, i);
            });
        }
        else if (_probDist == UNIFORM)
        {
            buildWRRBatch(numSamples, prevSize, [this](uint32_t u, size_t i) {
                BuildOneRRsetConstant(u, i);
            });
        }
        else if (_probDist == SKEWED || _probDist == WEIGHTS)
        {
            buildWRRBatch(numSamples, prevSize, [this](uint32_t u, size_t i) {
                BuildOneRRsetSkewed(u, i);
            });
        }
        else
        {
            buildWRRBatch(numSamples, prevSize, [this](uint32_t u, size_t i) {
                BuildOneRRset(u, i);
            });
        }

        _wrrBuildActive = false;
    }

    double EvalHyperedgeAvg()
    {
        return _hyperedgeAvgEval;
    }


    void display_hyperedge_stat()
    {
        size_t total_hyperedges = 0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            total_hyperedges += _RRsets[i].size();
        }

        double ave_hyperedge_size = 1.0 * total_hyperedges / _numRRsets;
        double var = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            double diff = _RRsets[i].size() - ave_hyperedge_size;
            var += (diff * diff);
        }

        var = var / _numRRsets;
        std::cout << "final average RR set size: " << ave_hyperedge_size << ", final variance: " << var << std::endl;
        return ;
    }

    double HyperedgeAvg()
    {
        size_t totalHyperedges = 0;
        double avgSize = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            totalHyperedges += _RRsets[i].size();
        }

        avgSize = 1.0 * totalHyperedges / _numRRsets;
        return avgSize;
    }

    double HyperedgeMedian()
    {
        std::vector<int> RRsetSize(_numRRsets);
        for (size_t i = 0; i < _numRRsets; i++)
        {
          RRsetSize[i]  = _RRsets[i].size();
          std::cout << RRsetSize[i] << std::endl;
        }

        std::sort(RRsetSize.begin(), RRsetSize.end());
        int index = _numRRsets/2;
        if (_numRRsets%2 == 1)
        {
            return RRsetSize[index];
        }
        else
        {
            return (RRsetSize[index] + RRsetSize[index-1])/2.0;
        }
    } 

    // ===== FRIM RR-set tracking (dummy node index = _numV) =====
    void FrimClearRRTracking()
    {
        for (auto i = _numV; i--;)
            FRset().swap(_FRsets[i]);
        FRset().swap(_FRsets[_numV]);

        _frimNumRR = 0;
        _frimRRweights.clear();
    }

    void FrimRegisterRRSample(const std::vector<uint32_t>& rr_nodes, double weight)
    {
        const size_t hyperIdx = _frimNumRR++;
        _frimRRweights.push_back(weight);

        for (uint32_t x : rr_nodes)
        {
            if (x < _numV)
                _FRsets[x].push_back(hyperIdx);
        }

        if (weight > 0.0)
            _FRsets[_numV].push_back(hyperIdx);
    }

    double FrimEstimateJFromDummy(double root_weight_sum, size_t num_rr) const
    {
        if (num_rr == 0 || root_weight_sum <= 0.0)
            return 0.0;

        double total = 0.0;
        for (size_t idx : _FRsets[_numV])
            total += _frimRRweights[idx];
        return (root_weight_sum / static_cast<double>(num_rr)) * total;
    }

    /// Reverse BFS RR-set: root skips xi; other nodes draw xi once per sample (cached);
    /// xi-fail nodes are never enqueued; q hit stops.
    FrimRRSample BuildOneFrimRRSample(
        const uint32_t v_root,
        const std::vector<double>& q,
        const std::vector<double>& xi) const
    {
        FrimRRSample sample;
        sample.root = v_root;
        sample.hit = false;
        sample.rr_nodes.clear();

        std::vector<uint8_t> visited(_numV, 0);
        // -1 unknown, 0 fail, 1 pass (one xi coin per node per RR-set)
        std::vector<int8_t> xi_pass(_numV, -1);
        std::queue<uint32_t> bfs_q;

        auto resolveXiPass = [&](uint32_t node) -> bool
        {
            if (node == v_root)
                return true;
            if (node >= xi.size())
                return false;
            if (xi_pass[node] >= 0)
                return xi_pass[node] == 1;
            const bool pass = dsfmt_gv_genrand_close_open() <= xi[node];
            xi_pass[node] = pass ? 1 : 0;
            return pass;
        };

        visited[v_root] = 1;
        bfs_q.push(v_root);

        while (!bfs_q.empty())
        {
            const uint32_t u = bfs_q.front();
            bfs_q.pop();

            if (u != v_root && !resolveXiPass(u))
                continue;

            sample.rr_nodes.push_back(u);

            if (u < q.size() && dsfmt_gv_genrand_open_close() <= q[u])
            {
                sample.hit = true;
                break;
            }

            for (const auto& nbr : _graph[u])
            {
                const uint32_t pred = nbr.first;
                if (visited[pred])
                    continue;
                if (dsfmt_gv_genrand_open_close() > nbr.second)
                    continue;
                if (pred != v_root && !resolveXiPass(pred))
                    continue;

                visited[pred] = 1;
                bfs_q.push(pred);
            }
        }

        return sample;
    }

    /// Reverse BFS without xi: save tree + per-node xi_gate for CRN re-evaluation.
    FrimRRStructureSample BuildOneFrimRRStructure(
        const uint32_t v_root,
        const std::vector<double>& q) const
    {
        FrimRRStructureSample sample;
        sample.root = v_root;
        sample.hit = false;
        sample.hit_node = 0;
        sample.parent.assign(_numV, UINT32_MAX);
        sample.xi_gate.assign(_numV, -1.0f);

        std::vector<uint8_t> visited(_numV, 0);
        std::queue<uint32_t> bfs_q;

        visited[v_root] = 1;
        sample.xi_gate[v_root] = 0.0f;
        sample.bfs_order.push_back(v_root);
        bfs_q.push(v_root);

        while (!bfs_q.empty())
        {
            const uint32_t u = bfs_q.front();
            bfs_q.pop();

            if (u < q.size() && dsfmt_gv_genrand_open_close() <= q[u])
            {
                sample.hit = true;
                sample.hit_node = u;
                break;
            }

            for (const auto& nbr : _graph[u])
            {
                const uint32_t pred = nbr.first;
                if (visited[pred])
                    continue;
                if (dsfmt_gv_genrand_open_close() > nbr.second)
                    continue;

                visited[pred] = 1;
                sample.parent[pred] = u;
                sample.xi_gate[pred] = static_cast<float>(dsfmt_gv_genrand_close_open());
                sample.bfs_order.push_back(pred);
                bfs_q.push(pred);
            }
        }

        return sample;
    }

    // Generate one RR set
    void BuildOneRRset(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
        const double gain = _wrrBuildActive ? wrrRootWeight(uStart) : (1.0 / (_vecOutdeg[uStart] + 1));
        for (auto node : _RRsets.back()) {
            _nodeGain[node] += gain;
        }
    }

    // Independent cascade with weighted probability
    void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            //if (_cascadeModel == IC)
            {
                if (_graph[expand].size() == 0) continue;

                double p =  _graph[expand][0].second;
                double log2Prob = Logarithm(1 - p);

                if (p < 1)
                {
                    double prob = dsfmt_gv_genrand_open_close();
                    int startPos = Logarithm(prob) / log2Prob;
                    int endPos = _graph[expand].size();

                    while (startPos < endPos)
                    {
                        const auto nbrId = _graph[expand][startPos].first;

                        if (_vecVisitBool[nbrId])
                        {
                            int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                            startPos += (increment + 1);
                            continue;
                        }

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                        int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                        startPos += increment + 1;
                    }
                }
                else
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
        const double gain = _wrrBuildActive ? wrrRootWeight(uStart) : (1.0 / _vecOutdeg[uStart]);
        for (auto& node : _RRsets.back()) {
            _nodeGain[node] += gain;
        }
    }

    /* independent cascade with constant probability */
    void BuildOneRRsetConstant(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        const double p =  _graph[0][0].second;
        const double const_prob = Logarithm(1 - p);

        if (p == 1)
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                //std::cout<<_graph[expand].size()<<std::endl;
                if (_graph[expand].size() == 0) continue;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    //std::cout<<nbr.first<<" "<<nbr.second<<std::endl;
                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }
        else
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    //std::cout<<"enter loop"<<std::endl;
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
        const double gain = _wrrBuildActive ? wrrRootWeight(uStart) : (1.0 / _vecOutdeg[uStart]);
        for (auto& node : _RRsets.back()) {
            _nodeGain[node] += gain;
        }
    }

    // independent cascade with skewed distribution
    void BuildOneRRsetSkewed(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        double p_threshold = 0.1;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            size_t out_degree = _graph[expand].size();

            if (out_degree > 0)
            {
                size_t startMin = 0;
                size_t endMax = out_degree;

                while (startMin < endMax)
                {
                    const auto &currentedge = _graph[expand][startMin];
                    const auto node_prob = currentedge.second;
                    const auto nbrId = currentedge.first;

                    if (node_prob < p_threshold)
                    {
                        break;
                    }

                    const auto randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > node_prob) continue;

                    if (_vecVisitBool[nbrId]) continue;

                    _vecVisitNode[numVisitNode++] =  nbrId;
                    _vecVisitBool[nbrId]  = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                while (startMin < endMax)
                {
                    double bucket_probability = _graph[expand][startMin].second;
                    const double log_prob = Logarithm(1 - bucket_probability);
                    double prob = dsfmt_gv_genrand_open_close();
                    startMin += floor(Logarithm(prob) / log_prob);

                    if (startMin >= endMax)
                    {
                        break;
                    }

                    const auto &currentedge = _graph[expand][startMin];
                    const auto nbrId = currentedge.first;
                    const auto accept_probability = currentedge.second;
                    double randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > accept_probability / bucket_probability || _vecVisitBool[nbrId])
                    {
                        continue;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
        const double gain = _wrrBuildActive ? wrrRootWeight(uStart) : (1.0 / _vecOutdeg[uStart]);
        for (auto& node : _RRsets.back()) {
            _nodeGain[node] += gain;
        }
    }

    // Evaluate the influence spread of a seed set on current generated RR sets
    double CalculateInf(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        for (auto seed : vecSeed)
        {
            for (auto node : _FRsets[seed])
            {
                vecBoolVst[node] = true;
            }
        }

        size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
        return 1.0 * count * _numV / _numRRsets;
    }

    double CalculateInf_W(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        double count = 0.0;

        for (auto seed : vecSeed)
        {
            for (auto RRidx : _FRsets[seed])
            {
                if (!vecBoolVst[RRidx]) {
                    vecBoolVst[RRidx] = true;
                    count += _RRweights[RRidx];
                }
            }
        }

        return count * _wrrNormConst / _numRRsets;
    }

    // Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
    double EfficInfVldtAlg(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
    {
        const double c = 2.0 * (exp(1.0) - 2.0);
        const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
        size_t numHyperEdge = 0;
        size_t numCoverd = 0;
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        while (numCoverd < LambdaL)
        {
            numHyperEdge++;
            size_t numVisitNode = 0, currIdx = 0;
            const auto uStart = dsfmt_gv_genrand_uint32_range(_numV);

            if (vecBoolSeed[uStart])
            {
                // Stop, this sample is covered
                numCoverd++;
                continue;
            }

            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    if (vecBoolSeed[nbrId])
                    {
                        // Stop, this sample is covered
                        numCoverd++;
                        goto postProcess;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                }
            }

postProcess:

            for (auto i = 0; i < numVisitNode; i++)
                _vecVisitBool[_vecVisitNode[i]] = false;
        }

        return 1.0 * numCoverd * _numV / numHyperEdge;
    }

    // Refresh the hypergraph
    void RefreshHypergraph()
    {
        if (_RRsets.size() != 0)
        {
            for (auto i = _numRRsets; i--;)
            {
                RRset().swap(_RRsets[i]);
            }

            RRsets().swap(_RRsets);

            for (auto i = _numV; i--;)
                FRset().swap(_FRsets[i]);
            FRset().swap(_FRsets[_numV]);
        }

        _numRRsets = 0;
        _hit = 0;
        _frimNumRR = 0;
        _frimRRweights.clear();
    }

    // Release memory
    void ReleaseMemory()
    {
        RefreshHypergraph();
        std::vector<bool>().swap(_vecVisitBool);
        Nodelist().swap(_vecVisitNode);
        FRsets().swap(_FRsets);
    }

    void BuildOneRRsetEarlyStopByVanilla(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            if (0 == _graph[expand].size())
            {
                continue;
            }

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }


    void BuildOneRRsetEarlyStopBySubsim(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            const double p =  _graph[expand][0].second;

            if (0 == _graph[expand].size())
            {
                continue;
            }

            if (p >= 1.0)
            {
                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                continue;
            }

            const double const_prob = Logarithm(1 - p);
            int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
            int endPos = _graph[expand].size();

            while (startPos < endPos)
            {
                const auto nbrId = _graph[expand][startPos].first;

                if (!_vecVisitBool[nbrId])
                {
                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }

                int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                startPos += increment + 1;
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    void BuildRRsetsEarlyStop(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR sets with early stop by Vanilla method" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopByVanilla(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
        else
        {
            // std::cout << "Sample RR sets with early stop By SUBSIM" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopBySubsim(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
    }

    double EvalSeedSetInfByVanilla(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInfBySubsim(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                if (p >= 1.0)
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;

                        if (connSet.find(nbrId) != connSet.end())
                        {
                            numCovered++;
                            goto finished;
                        }
                    }

                    continue;
                }

                const double const_prob = Logarithm(1 - p);
                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                    }

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInf(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        if (_isVanilla)
        {
            return EvalSeedSetInfByVanilla(connSet, numSamples);
        }
        else
        {
            return EvalSeedSetInfBySubsim(connSet, numSamples);
        }
    }

    double CalculateInfEarlyStop()
    {
        return 1.0 * _hit * _numV / _numRRsets;
    }
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
