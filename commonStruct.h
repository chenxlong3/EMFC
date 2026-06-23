#pragma once
/// Node list
typedef std::vector<uint32_t> Nodelist;
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

/// CRN RR structure: reverse BFS tree without xi filtering; xi_gate fixed per node.
struct FrimRRStructureSample
{
    uint32_t root = 0;
    bool hit = false;
    uint32_t hit_node = 0;
    std::vector<uint32_t> bfs_order;
    std::vector<uint32_t> parent;
    std::vector<float> xi_gate;
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
enum FuncType {FORMAT, IM, FUNC_ERROR, WIM, SUBSIM, OUTDEG, PROB, RAND, AIS, ScaLIM, ADIL, GREEDY, CHANGE_PROB, DG, EVAL, PREP_CAND, STAT, FRIM_RR, FRIM_RR_NAIVE, FRIM_RR_CRN, FRIM_MC, FRIM_MC_NAIVE};
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
