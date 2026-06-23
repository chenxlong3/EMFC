#pragma once

struct FrimLiveEdgeSample
{
    std::vector<double> seed_gates;
    std::vector<std::vector<uint32_t>> kept_out;
    std::vector<double> gates;
    std::vector<uint8_t> active;
};

struct FrimXiResult
{
    std::vector<double> xi;
    std::vector<double> alpha;
    double J_hat = 0.0;
    double J_method_rr = 0.0;
    double J_method_mc = 0.0;
    double J_method_mc_naive = 0.0;
    double J_method_rr_naive = 0.0;
    double J_method_rr_crn = 0.0;
    FrimRunInfo run_info;
};

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

    // ===== FRIM xi selection (implement-spec.md) =====
    static std::vector<double> buildCumulativeWeights(const std::vector<double>& weights);
    static uint32_t sampleByCumulativeWeights(const std::vector<double>& cumWeights);
    static double nodeObjectiveValue(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& xi,
        uint32_t nodeId);

    FrimRRSample frimSampleOneRR(
        const std::vector<double>& cum_root_weights,
        const std::vector<double>& q,
        const std::vector<double>& xi) const;
    FrimRRSample frimReverseBfsFirstHitRs(
        uint32_t v_root,
        const std::vector<double>& q,
        const std::vector<double>& xi) const;
    std::vector<FrimRRSample> frimBuildRRSamples(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi) const;
    double frimEstimateJ(
        const std::vector<double>& xi,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        size_t num_rr);
    void frimMethodRR(
        const std::vector<FrimRRSample>& samples,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double root_weight_sum,
        size_t num_rr,
        double xi_lo,
        std::vector<double>& xi,
        std::vector<double>& alpha) const;
    void frimComputeActive(
        const FrimLiveEdgeSample& sample,
        const std::vector<double>& q,
        const std::vector<double>& xi,
        std::vector<uint8_t>& active,
        uint32_t silent_node = UINT32_MAX) const;
    std::vector<uint8_t> frimForwardReachFrom(
        const FrimLiveEdgeSample& sample,
        const std::vector<double>& xi,
        uint32_t src) const;
    double frimObjectiveOnSample(
        const FrimLiveEdgeSample& sample,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& xi) const;
    double frimFlipDelta(
        const FrimLiveEdgeSample& sample,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& xi,
        uint32_t u,
        double xi_new) const;
    void frimApplyFlip(
        FrimLiveEdgeSample& sample,
        const std::vector<double>& q,
        const std::vector<double>& xi,
        uint32_t u,
        double xi_new,
        double xi_lo) const;
    FrimXiResult frimMethodMC(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_mc,
        int max_sweeps,
        double delta_tol) const;
    FrimLiveEdgeSample frimSampleLiveEdge(
        const Graph& forwardGraph) const;
    double frimObjectiveWithXi(
        const FrimLiveEdgeSample& sample,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& xi) const;
    FrimXiResult frimMethodMCNaive(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_mc,
        int max_sweeps,
        double delta_tol) const;
    FrimXiResult frimMethodRRNaive(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_rr,
        int max_sweeps,
        double delta_tol);
    FrimRRStructureSample frimBuildOneRRStructure(
        const std::vector<double>& cum_tau,
        const std::vector<double>& q) const;
    std::vector<FrimRRStructureSample> frimBuildRRStructureSamples(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q) const;
    static double frimRRStructureSampleWeight(
        const FrimRRStructureSample& sample,
        const std::vector<double>& xi,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        uint32_t override_node = UINT32_MAX,
        double override_xi = 0.0);
    double frimEstimateJFromRRStructures(
        const std::vector<FrimRRStructureSample>& samples,
        const std::vector<double>& xi,
        const std::vector<double>& tau,
        const std::vector<double>& lam) const;
    FrimXiResult frimMethodRRCrn(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_rr,
        int max_sweeps,
        double delta_tol);
    // ===== end FRIM xi selection =====

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

    /// FRIM RR method: tau-weighted RR one-pass scan.
    FrimXiResult frim_solve_rr(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo = 0.5,
        size_t num_rr = 10000);

    /// Naive RR: each sweep re-estimates J with xi_u in {xi_lo, 1} (fresh R per estimate).
    FrimXiResult frim_solve_rr_naive(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_rr = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);

    /// CRN RR: tau-root structural samples reused; per-node {xi_lo,1} via fixed xi_gate coins.
    FrimXiResult frim_solve_rr_crn(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_rr = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);

    /// RR J_hat at fixed xi (resamples num_rr sets).
    double frim_estimate_j_at_xi(
        const std::vector<double>& xi,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        size_t num_rr);

    /// MC gate-objective at fixed xi (live-edge CRN-style average, num_mc samples).
    double frim_estimate_mc_gate_at_xi(
        const Graph& forwardGraph,
        const std::vector<double>& xi,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        size_t num_mc) const;

    /// FRIM MC method: CRN live-edge coordinate ascent.
    FrimXiResult frim_solve_mc(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_mc = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);

    /// Naive MC: resample live-edge graphs each sweep; full BFS per flip test.
    FrimXiResult frim_solve_mc_naive(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_mc = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;
