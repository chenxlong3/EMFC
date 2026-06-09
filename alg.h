#pragma once

class Alg
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    /// Upper bound in the last round for __mode=1.
    double _boundLast = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2.
    double _boundMin = DBL_MAX;
    /// Two hyper-graphs, one is used for selecting seeds and the other is used for validating influence.
    THyperGraph _hyperGraph, _hyperGraphVldt;
    /// Result object.
    TResult& _res;
    /// Seed set.
    Nodelist _vecSeed;
    ProbDist _probDist = WC;
    size_t _norm_const = 0;

    double _baseNumRRsets = 0.0;

    std::vector<uint32_t> _vecOutDegree;
    std::vector<uint32_t> _vecVldtInf;

    /// Maximum coverage by lazy updating.
    double MaxCoverVanilla(const int targetSize);
    double MaxCoverOutDegPrority(const int targetSize);
    double MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize);
    double MaxCoverSentinelSet(const int targetSize, const int totalTargetSize);

    /// Maximum coverage by maintaining the top-k marginal coverage.
    double MaxCoverTopK(const int targetSize);
    /// Maximum coverage.
    double MaxCover(const int targetSize);
    double MaxWCover(const int targetSize);


public:
    Alg(const Graph& graph, TResult& tRes) : _hyperGraph(graph), _hyperGraphVldt(graph), _res(tRes)
    {
        _numV = _hyperGraph.get_nodes();
        _numE = _hyperGraph.get_edges();
        _norm_const = _numE;
        _vecOutDegree = std::vector<uint32_t>(_numV);

        for (auto &nbrs : graph)
        {
            for (auto &node : nbrs)
            {
                _vecOutDegree[node.first]++;
            }
        }
    }

    ~Alg()
    {
    }

    /// Set cascade model.
    void set_prob_dist(const ProbDist weight);
    void set_vanilla_sample(const bool isVanilla);
    void set_wrr_sample_mode(const WRRSampleMode mode);

    void RefreshHypergraph()
    {
        _hyperGraph.RefreshHypergraph();
        _hyperGraphVldt.RefreshHypergraph();
    }
    /// Evaluate influence spread for the seed set constructed
    double EfficInfVldtAlg();
    /// Evaluate influence spread for a given seed set
    double EfficInfVldtAlg(const Nodelist vecSeed);


    double estimateRRSize();

    double subsimOnly(const int targetSize, const double epsilon, const double delta);
    double subsimWeight(const int targetSize, const double epsilon, const double delta);
    double subsimWithTrunc(const int targetSize, const double epsilon, const double delta);
    double IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr);

    double FindFixSub(const int targetSize, const int totalTargetSize, const double epsilon, const double delta);
    double FindRemSet(const int targetSize, const double epsilon, const double targeEpsilon, const double delta);
    double FindDynamSub(const int totalTargetSize, const double epsilon, const double delta);

    double subsimWithHIST(const int targetSize, const double epsilon, const double delta);
    
    // Fixed number of RR sets sampling and then run greedy algorithm
    double fixed_subsim(const int targetSize, int num_samples);
    double fixed_subsimW(const int targetSize, int num_samples);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;
