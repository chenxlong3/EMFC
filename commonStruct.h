#pragma once
#include <unordered_map>
/// Node list
typedef std::vector<uint32_t> Nodelist;
/// Sparse forward adjacency: only nodes with out-edges appear as keys.
using FrimSparseFwdAdj = std::unordered_map<uint32_t, std::vector<uint32_t>>;
/// Sparse per-node byte flags (xi_active, reachability, Tarjan marks).
using FrimByteMap = std::unordered_map<uint32_t, uint8_t>;
/// Sparse per-node int state (Tarjan disc/low/parent).
using FrimIntMap = std::unordered_map<uint32_t, int>;
/// CRN xi coins in [0,1] for nodes in an RR-graph sample.
using FrimXiGateMap = std::unordered_map<uint32_t, float>;
/// Edge structure, neighbor id and the edge weight
typedef std::pair<uint32_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
typedef std::vector<size_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
double decay_factor = 1.0;

/// Per-node hyperparameters loaded from {graph}.nodehyper.vec at runtime.
struct NodeHyperParamsData
{
    std::vector<double> q;
    std::vector<double> tau;
    std::vector<double> lam;
    bool loaded = false;
};
NodeHyperParamsData g_nodeHyperParams;

/// Loaded FRIM xi vector for benefit evaluation.
struct FrimXiData
{
    std::vector<double> xi;
    double J_hat = 0.0;
    double J_method_rr = 0.0;
    double J_method_mc = 0.0;
    bool loaded = false;
};
FrimXiData g_frimXi;

/// One FRIM RR sample from reverse BFS with live edges, xi forwarding, and q seed hit.
struct FrimRRSample
{
    uint32_t root = 0;
    bool hit = false;
    std::vector<uint32_t> rr_nodes;
};

/// RR sample with full forward edge list on the sampled RR-graph (for K / Tarjan).
struct FrimRRGraphSample
{
    uint32_t root = 0;
    /// Virtual seed id (= numV) used only for Tarjan; not stored in nodes/adj.
    uint32_t platform_id = UINT32_MAX;
    bool hit = false;
    /// q-hit seed nodes recorded during reverse BFS.
    std::vector<uint32_t> hit_nodes;
    std::vector<uint32_t> nodes;
    /// Forward live edges pred -> u within the RR set.
    FrimSparseFwdAdj forward_adj;
    /// CRN xi coins in [0,1] for nodes in this RR-graph (sparse; root omitted).
    FrimXiGateMap xi_gate;
};

/// CRN RR structure: full reverse-BFS subgraph + fixed xi_gate per node.
struct FrimRRStructureSample
{
    uint32_t root = 0;
    bool hit = false;
    /// First q-hit node (dequeue order); kept for logging / legacy.
    uint32_t hit_node = 0;
    /// All q-hit nodes; q-hit stops expansion at u only, BFS continues elsewhere.
    std::vector<uint32_t> hit_nodes;
    /// All nodes in the sampled RR subgraph.
    std::vector<uint32_t> nodes;
    /// Live forward edges (pred -> u) within the RR subgraph.
    FrimSparseFwdAdj forward_adj;
    /// CRN xi coins in [0,1] for nodes in this subgraph (sparse; root omitted or 0).
    FrimXiGateMap xi_gate;
};

/// Runtime stats for a FRIM solve run (written to result/info).
struct FrimRunInfo
{
    std::string method;
    int rand_seed = 0;
    size_t num_v = 0;
    size_t num_rr = 0;
    size_t num_mc = 0;
    int max_sweeps = 0;
    int sweeps_completed = 0;
    double xi_lo = 0.5;
    size_t num_xi_one = 0;
    size_t num_xi_lo = 0;
    double time_total_sec = 0.0;
    double time_sample_sec = 0.0;
    double time_solve_sec = 0.0;
    double time_estimate_sec = 0.0;
};
/// One reverse reachable set
typedef std::vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRsets;
bool optflag;
enum ProbDist {WEIGHTS, UNIFORM, WC, WC2, SKEWED, PROB_DIST_ERROR};
enum FuncType {FORMAT, HYP_DERIVE, IM, FUNC_ERROR, WIM, SUBSIM, OUTDEG, PROB, RAND, AIS, ScaLIM, ADIL, GREEDY, CHANGE_PROB, DG, EVAL, PREP_CAND, STAT, FRIM_RR, FRIM_RR_NAIVE, FRIM_RR_GRAPH, FRIM_RR_GRAPH_PRUNE, FRIM_RR_GRAPH_ROOT_STAT, FRIM_SUBSIM, FRIM_PRUNE, FRIM_PRUNE_LO, FRIM_PRUNE_BOTH, FRIM_MC_CRN, FRIM_MC_NAIVE, FRIM_OUTNAME};
/// Weighted RR root sampling mode for subsimW (-func=wim)
enum WRRSampleMode { WRR_INDEG, WRR_WINDEG, WRR_WOUTDEG, WRR_OUTDEG, WRR_UNIFORM, WRR_SAMPLE_ERROR };
enum CascadeModel { IC, LT };

/// Node element with id and a property value
typedef struct NodeElement
{
    int id;
    double value;
} NodeEleType;

/// Smaller operation for node element
struct smaller
{
    bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
    {
        return (Ele1.value < Ele2.value);
    }
};

/// Greater operation for node element
struct greater
{
    bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
    {
        return (Ele1.value > Ele2.value);
    }
};
