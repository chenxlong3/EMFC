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
    double J_method_rr_graph = 0.0;
    FrimRunInfo run_info;
};

/// K-L gradient prune analysis at xi ≡ 1 (L_ub = lam_u * q-hit count at xi≡1, root u).
struct FrimPruneResult
{
    std::vector<double> K;
    std::vector<double> K_lb;
    std::vector<double> L_ub;
    std::vector<double> alpha_lb;
    std::vector<uint8_t> pruned_to_one;
    size_t num_pruned = 0;
    size_t num_uncertain = 0;
    size_t num_rr = 0;
    double tau_sum = 0.0;
    double hoeffding_margin = 0.0;
};

/// K upper / L lower prune analysis at xi ≡ xi_lo (L_lb = lam_u * q-hit count at xi≡xi_lo, root u).
struct FrimPruneLoResult
{
    std::vector<double> K;
    std::vector<double> K_ub;
    std::vector<double> L_lb;
    std::vector<double> alpha_ub;
    std::vector<uint8_t> pruned_to_lo;
    size_t num_pruned = 0;
    size_t num_uncertain = 0;
    size_t num_rr = 0;
    double tau_sum = 0.0;
    double hoeffding_margin = 0.0;
    double xi_lo = 0.5;
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
    double frimEstimateBenefitJ(
        const Graph& forwardGraph,
        const std::vector<double>& xi,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        CascadeModel model,
        size_t num_mc) const;
    FrimXiResult frimMethodMCNaive(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_mc,
        int max_sweeps,
        double delta_tol,
        CascadeModel model) const;
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
    std::vector<FrimRRSample> frimBuildTauRootRRSamplesAtXiOne(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q) const;
    std::vector<FrimRRGraphSample> frimBuildTauRootRRGraphSamplesAtXiOne(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q) const;
    std::vector<FrimRRGraphSample> frimBuildTauRootRRGraphSamples(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q,
        const std::vector<double>& xi) const;
    std::vector<FrimRRGraphSample> frimBuildTauRootRRGraphSamplesCrn(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q,
        bool store_hit_only = true,
        size_t* num_discarded_out = nullptr) const;
    static bool frimRRGraphNodeBlocks(
        const FrimRRGraphSample& sample,
        uint32_t u);
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
    FrimXiResult frimMethodRRGraph(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo,
        size_t num_rr,
        int max_sweeps,
        double delta_tol,
        bool combined_prune_init = false,
        double eps = 0.1,
        double hoeffding_delta = 0.05,
        double hoeffding_margin_scale = 1.0,
        double hoeffding_margin_override = -1.0,
        bool rr_graph_gate_sweep_index = true,
        bool store_hit_only = true,
        const std::vector<uint32_t>* subsim_fix_one = nullptr);
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

    /// Naive RR: J_stay vs J_flip (one frimEstimateJ each); J_hat = J_flip on accept.
    FrimXiResult frim_solve_rr_naive(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_rr = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);

    /// RR-graph: tau-root structural samples reused; per-node {xi_lo,1} via fixed xi_gate coins.
    FrimXiResult frim_solve_rr_graph(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_rr = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9,
        bool rr_graph_gate_sweep_index = true,
        bool store_hit_only = true);

    /// Count tau-root assignments over R RR-graph sampling attempts; print min/max stats.
    void frimStatAndPrintRRGraphRootCounts(
        size_t num_rr,
        const std::vector<double>& tau,
        const std::vector<double>& q,
        bool store_hit_only,
        bool build_rr_graphs) const;

    /// SubSim heuristic: top-(n/100) IM seeds get xi=1, rest xi=xi_lo; no FRIM estimate.
    FrimXiResult frim_solve_subsim(
        double xi_lo = 0.5,
        double subsim_eps = 0.5,
        double subsim_delta = 0.0);

    /// RR-graph + combined K/L prune init; sweep only non-fixed nodes.
    FrimXiResult frim_solve_rr_graph_prune(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_rr = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9,
        double eps = 0.1,
        double hoeffding_delta = 0.05,
        double hoeffding_margin_scale = 1.0,
        double hoeffding_margin_override = -1.0,
        bool rr_graph_gate_sweep_index = true,
        bool store_hit_only = true);

    /// At xi≡1: K_lb from xi_lo pass-through/block samples; L_ub=lam*hit count; prune if K_lb-L_ub>0.
    FrimPruneResult frimAnalyzeKMinusLPrune(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        size_t num_rr,
        double eps = 0.1,
        double hoeffding_delta = 0.05) const;

    FrimPruneResult frimAnalyzeKMinusLPruneFromSamples(
        const std::vector<FrimRRGraphSample>& samples,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        double eps = 0.1,
        double hoeffding_delta = 0.05,
        double hoeffding_margin_scale = 1.0,
        double hoeffding_margin_override = -1.0) const;

    void frimRunKMinusLPruneAnalysis(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        size_t num_rr,
        double eps = 0.1,
        double hoeffding_delta = 0.05) const;

    /// At xi≡xi_lo: K_ub excludes pass-through samples; L_lb=lam*hit count at xi≡xi_lo; prune if K_ub-L_lb<0.
    FrimPruneLoResult frimAnalyzeKUbLlbAtXiLo(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        size_t num_rr,
        double eps = 0.1,
        double hoeffding_delta = 0.05) const;

    FrimPruneLoResult frimAnalyzeKUbLlbAtXiLoFromSamples(
        const std::vector<FrimRRGraphSample>& samples,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        double eps = 0.1,
        double hoeffding_delta = 0.05,
        double hoeffding_margin_scale = 1.0,
        double hoeffding_margin_override = -1.0) const;

    void frimRunKUbLlbPruneAnalysis(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        size_t num_rr,
        double eps = 0.1,
        double hoeffding_delta = 0.05) const;

    void frimRunCombinedPruneAnalysis(
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        double xi_lo,
        size_t num_rr,
        int rand_seed,
        double eps = 0.1,
        double hoeffding_delta = 0.05) const;

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

    /// FRIM MC-CRN: CRN live-edge coordinate ascent.
    FrimXiResult frim_solve_mc_crn(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_mc = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9);

    /// Naive MC: one benefit_inf_eval each for J_stay and J_flip; J_hat = J_flip on accept.
    FrimXiResult frim_solve_mc_naive(
        const Graph& forwardGraph,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const std::vector<double>& q,
        const std::vector<double>& xi_init,
        double xi_lo = 0.5,
        size_t num_mc = 1000,
        int max_sweeps = 10,
        double delta_tol = 1e-9,
        CascadeModel model = IC);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;
