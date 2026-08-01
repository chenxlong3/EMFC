#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

void init_random_seed(int rand_seed)
{
    const uint32_t seed = (rand_seed != 0)
        ? static_cast<uint32_t>(rand_seed)
        : static_cast<uint32_t>(time(nullptr));
    dsfmt_gv_init_gen_rand(seed);
    if (rand_seed != 0)
        std::cout << "Random seed: " << rand_seed << std::endl;
}

static size_t frimResolveRRGraphNumRR(TArgument& Arg, size_t num_v, double tau_sum, const char* tag)
{
    const size_t num_rr = Arg.resolveFrimRRGraphNumSamples(num_v, tau_sum);
    if (!Arg._frim_rr_fixed)
    {
        const double n = static_cast<double>(num_v);
        const double base_delta = (Arg._delta > 0.0) ? Arg._delta : (1.0 / n);
        const double delta = base_delta / 2.0;
        std::cout << tag
                  << " num_rr = (8+2*eps)*T*(log(1/delta)+n*log(2)+log(2))/(eps^2*n): n="
                  << num_v
                  << " T=" << tau_sum
                  << " delta=" << delta
                  << " eps=" << Arg._eps
                  << " -> R=" << num_rr << std::endl;
    }
    Arg._frim_rr = num_rr;
    return num_rr;
}

int main(int argc, char* argv[])
{

    TArgument Arg(argc, argv);

    if (Arg._probDist == PROB_DIST_ERROR)
    {
        LogInfo("The input probability distribution is not supported:", Arg._probDistStr);
        LogInfo("The supported probability distribution: weights, wc, uniform, skewed");
        return 0;
    }

    if (Arg._func == FUNC_ERROR)
    {
        LogInfo("The input func is not supported: ", Arg._funcStr);
        LogInfo("The supported func: format, hyp_derive, im, wim, SUBSIM, OUTDEG, PROB, RAND, AIS, ScaLIM, ADIL, GREEDY, DG, EVAL, PREP_CAND, STAT, frim_rr, frim_rr_naive, frim_rr_graph, frim_rr_graph_prune, frim_rr_root_stat, frim_subsim, frim_prune, frim_prune_lo, frim_prune_both, frim_mc_crn, frim_mc_naive, frim_outname");
        return 0;
    }

    if (Arg._func == FRIM_OUTNAME)
    {
        const FuncType method = Arg.decode_frim_eval_method();
        if (method == FUNC_ERROR)
        {
            std::cout << "Error: frim_outname requires -method=frim_rr, frim_rr_naive, "
                      << "frim_rr_graph, frim_rr_graph_prune, frim_mc_crn, or frim_mc_naive"
                      << std::endl;
            return 1;
        }
        Arg.build_frim_outfilename(method);
        std::cout << Arg._outFileName << std::endl;
        return 0;
    }

    if (Arg._func == WIM && Arg._wrrSampleMode == WRR_SAMPLE_ERROR)
    {
        LogInfo("The input wrr_sample is not supported:", Arg._wrrSampleStr);
        LogInfo("The supported wrr_sample: indeg, windeg, woutdeg, outdeg, uniform");
        return 0;
    }

    if (Arg._func == FORMAT && Arg._lamModeError)
    {
        LogInfo("The input lam_mode is not supported:", Arg._lamModeStr);
        LogInfo("The supported lam_mode: uniform, unirand, two_tier, exponential");
        return 0;
    }

    if ((Arg._func == FORMAT || Arg._func == HYP_DERIVE) && Arg._qModeError)
    {
        LogInfo("The input q_mode is not supported:", Arg._qModeStr);
        LogInfo("The supported q_mode: active_inactive, exponential, unirand");
        return 0;
    }

    if (Arg._func == FORMAT && Arg._tauModeError)
    {
        LogInfo("The input tau_mode is not supported:", Arg._tauModeStr);
        LogInfo("The supported tau_mode: exponential, ranked, soft, q_normal, uniform, exp_random");
        return 0;
    }

    if (Arg._func != HYP_DERIVE)
        init_random_seed(Arg._rand_seed);

    const std::string infilename = Arg._dir + "/" + Arg._graphname;
    if (Arg._func == FORMAT)
    {
        GraphBase::FormatGraph(
            infilename,
            Arg._probDist,
            Arg._wcVar,
            Arg._probEdge,
            Arg._skewType,
            GraphBase::makeNodeHyperParamsConfig(
                Arg._lamModeStr,
                Arg._lamUniformValue,
                Arg._lamUnirandLo,
                Arg._lamUnirandHi,
                Arg._activeUserRatio,
                Arg._activeQMean,
                Arg._activeQVar,
                Arg._inactiveQMean,
                Arg._inactiveQVar,
                Arg._qModeStr,
                Arg._qExpScale,
                Arg._qUnirandLo,
                Arg._qUnirandHi,
                Arg._tauModeStr,
                Arg._tauLo,
                Arg._tauHi,
                Arg._tauJitter,
                Arg._tauQVar));
        return 0;
    }

    if (Arg._func == HYP_DERIVE)
    {
        GraphBase::HypDeriveSpec spec;
        if (!Arg._hypProfilePath.empty())
        {
            if (!GraphBase::loadHypDeriveProfile(Arg._hypProfilePath, spec))
                return 1;
        }
        else
        {
            if (Arg._hypOutputSuffix.empty() && Arg._nodehyperSuffix.empty())
            {
                std::cout << "Error: hyp_derive requires -hyp_output_suffix or -nodehyper_suffix"
                          << " (or -hyp_profile)" << std::endl;
                return 1;
            }
            spec.base_suffix = Arg._nodehyperBaseSuffix;
            spec.output_suffix = !Arg._hypOutputSuffix.empty()
                ? Arg._hypOutputSuffix
                : Arg._nodehyperSuffix;
            spec.keep_tau = true;
            spec.keep_lam = true;
            if (Arg._qModeStr != "keep")
            {
                spec.override_q = true;
                spec.q_config = GraphBase::makeNodeHyperParamsConfig(
                    Arg._lamModeStr,
                    Arg._lamUniformValue,
                    Arg._lamUnirandLo,
                    Arg._lamUnirandHi,
                    Arg._activeUserRatio,
                    Arg._activeQMean,
                    Arg._activeQVar,
                    Arg._inactiveQMean,
                    Arg._inactiveQVar,
                    Arg._qModeStr,
                    Arg._qExpScale,
                    Arg._qUnirandLo,
                    Arg._qUnirandHi,
                    Arg._tauModeStr,
                    Arg._tauLo,
                    Arg._tauHi,
                    Arg._tauJitter,
                    Arg._tauQVar);
            }
            spec.rand_seed = Arg._rand_seed;
        }

        const int derive_seed = (spec.rand_seed >= 0) ? spec.rand_seed : Arg._rand_seed;
        init_random_seed(derive_seed);

        std::cout << "HYP_DERIVE: base="
                  << (spec.base_suffix.empty() ? "default" : spec.base_suffix)
                  << " -> output=" << spec.output_suffix
                  << ", override_q=" << (spec.override_q ? "yes" : "no")
                  << ", keep_tau=" << (spec.keep_tau ? "yes" : "no")
                  << ", keep_lam=" << (spec.keep_lam ? "yes" : "no")
                  << std::endl;

        if (!GraphBase::deriveAndSaveNodeHyperParams(infilename, spec))
            return 1;
        return 0;
    }

    if (!TIO::LoadGraphNodeHyperParamsData(
            infilename, g_nodeHyperParams, Arg._nodehyperSuffix))
    {
        std::cout << "Warning: node hyper params file not found for " << infilename
                  << ", run -func=format first." << std::endl;
    }

    if (Arg._func == STAT)
    {
        Graph graph;
        GraphBase::LoadForwardGraph(graph, infilename);
        GraphBase::computeAndAppendDegreeStatsToAttr(graph, Arg._graphname, Arg._dir);
        std::cout << "STAT: appended degree stats to " << Arg._dir << "/" << Arg._graphname << ".attr" << std::endl;
        return 0;
    }

    // Common variables
    int seedSize = Arg._seedsize;

    Timer mainTimer("main");
    // Forward/reverse graphs are saved at format time (.vec.graph / .vec.rvs.graph).
    Graph graph;
    Graph forwardGraph;
    if (Arg._func == EVAL)
    {
        GraphBase::LoadForwardGraph(graph, infilename);
    }
    else if (Arg._func == FRIM_MC_CRN || Arg._func == FRIM_MC_NAIVE)
    {
        GraphBase::LoadGraph(graph, infilename);
        GraphBase::LoadForwardGraph(forwardGraph, infilename);
    }
    else
    {
        GraphBase::LoadGraph(graph, infilename);
    }
    int probDist = GraphBase::LoadGraphProbDist(infilename);

    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();

    // FRIM RR method: one-pass RR scan, save xi.
    if (Arg._func == FRIM_RR)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_RR);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_rr(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            Arg._xi_lo,
            Arg._frim_rr);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_rr";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_rr = Arg._frim_rr;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_naive,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM RR naive: resample RR sets per flip test, save xi.
    if (Arg._func == FRIM_RR_NAIVE)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        const double xi_start = Arg._xi_init_lo ? Arg._xi_lo : 1.0;
        std::vector<double> xi_init(graph.size(), xi_start);
        std::cout << "FRIM-RR-naive xi init: all " << xi_start << std::endl;

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_RR_NAIVE);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_rr_naive(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            xi_init,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._frim_max_sweeps);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_rr_naive";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_rr = Arg._frim_rr;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_naive,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM RR-graph: fixed tau-root structural samples, xi_gate CRN per node.
    if (Arg._func == FRIM_RR_GRAPH)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        const double xi_start = Arg._xi_init_lo ? Arg._xi_lo : 1.0;
        std::vector<double> xi_init(graph.size(), xi_start);
        std::cout << "FRIM-RR-GRAPH xi init: all " << xi_start << std::endl;
#if 0  // warm-start disabled
        if (Arg._rr_graph_warm_start)
        {
            std::cout << "FRIM-RR-GRAPH warm-start: first R/"
                      << Arg._rr_graph_warm_start_div
                      << " samples + " << Arg._rr_graph_warm_start_sweeps
                      << " sweep(s), then full R" << std::endl;
        }
#endif

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        double tau_sum = 0.0;
        for (double t : g_nodeHyperParams.tau) tau_sum += t;
        frimResolveRRGraphNumRR(Arg, graph.size(), tau_sum, "FRIM-RR-GRAPH");
        Arg.build_frim_outfilename(FRIM_RR_GRAPH);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_rr_graph(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            xi_init,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._frim_max_sweeps,
            1e-9,
            Arg._rr_graph_gate_sweep_index,
            Arg._rr_graph_store_hit_only);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_rr_graph";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_rr = Arg._frim_rr;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_graph,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // Count tau-root assignments over R RR-graph sampling attempts.
    if (Arg._func == FRIM_RR_GRAPH_ROOT_STAT)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        double tau_sum = 0.0;
        for (double t : g_nodeHyperParams.tau) tau_sum += t;
        frimResolveRRGraphNumRR(Arg, graph.size(), tau_sum, "FRIM-RR-ROOT-STAT");

        tAlg.frimStatAndPrintRRGraphRootCounts(
            Arg._frim_rr,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.q,
            Arg._rr_graph_store_hit_only,
            Arg._rr_root_stat_build);
        return 0;
    }

    // SubSim heuristic: IM seeds xi=1, rest xi=xi_lo (no FRIM estimate).
    if (Arg._func == FRIM_SUBSIM)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_SUBSIM);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        const double subsim_delta = (Arg._delta > 0.0)
            ? Arg._delta
            : (1.0 / static_cast<double>(graph.size()));

        FrimXiResult frimRes = tAlg.frim_solve_subsim(
            Arg._xi_lo,
            Arg._eps,
            subsim_delta);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_subsim";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_rr = 0;
        frimRes.run_info.max_sweeps = 0;
        frimRes.run_info.xi_lo = Arg._xi_lo;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_graph,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM RR-graph + combined K/L prune init, then sweep active nodes only.
    if (Arg._func == FRIM_RR_GRAPH_PRUNE)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_init(graph.size(), 1.0);
        std::cout << "FRIM-RR-GRAPH-PRUNE xi init: all 1.0, then K/L prune fixes nodes"
                  << std::endl;
#if 0  // warm-start disabled
        if (Arg._rr_graph_warm_start)
        {
            std::cout << "FRIM-RR-GRAPH-PRUNE warm-start: first R/"
                      << Arg._rr_graph_warm_start_div
                      << " samples + " << Arg._rr_graph_warm_start_sweeps
                      << " sweep(s), then full R + prune" << std::endl;
        }
#endif

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        double tau_sum = 0.0;
        for (double t : g_nodeHyperParams.tau) tau_sum += t;
        frimResolveRRGraphNumRR(Arg, graph.size(), tau_sum, "FRIM-RR-GRAPH-PRUNE");
        Arg.build_frim_outfilename(FRIM_RR_GRAPH_PRUNE);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_rr_graph_prune(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            xi_init,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._frim_max_sweeps,
            1e-9,
            Arg._eps,
            Arg._hoeffding_delta,
            Arg._hoeffding_margin_scale,
            Arg._hoeffding_margin,
            Arg._rr_graph_gate_sweep_index,
            Arg._rr_graph_store_hit_only);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_rr_graph_prune";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_rr = Arg._frim_rr;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_graph,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM K-L prune analysis at xi ≡ 1.
    if (Arg._func == FRIM_PRUNE)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        std::cout << "---The Begin of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        tAlg.frimRunKMinusLPruneAnalysis(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._eps);
        std::cout << "---The End of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        return 0;
    }

    // FRIM K_ub / L_lb prune analysis at xi ≡ xi_lo.
    if (Arg._func == FRIM_PRUNE_LO)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        std::cout << "---The Begin of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_lo_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        tAlg.frimRunKUbLlbPruneAnalysis(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._eps);
        std::cout << "---The End of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_lo_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        return 0;
    }

    // Combined K_lb/L_ub (xi=1) and K_ub/L_lb (xi=xi_lo) prune analysis.
    if (Arg._func == FRIM_PRUNE_BOTH)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        std::cout << "---The Begin of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_both_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        tAlg.frimRunCombinedPruneAnalysis(
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            Arg._xi_lo,
            Arg._frim_rr,
            Arg._rand_seed,
            Arg._eps);
        std::cout << "---The End of " << Arg._graphname << "_"
                  << Arg._rand_seed << "_frim_prune_both_" << Arg._probDistStr
                  << "_R" << Arg._frim_rr << "---\n";
        return 0;
    }

    // FRIM MC-CRN: CRN live-edge coordinate ascent, save xi.
    if (Arg._func == FRIM_MC_CRN)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        const double xi_start = Arg._xi_init_lo ? Arg._xi_lo : 1.0;
        std::vector<double> xi_init(graph.size(), xi_start);
        std::cout << "FRIM-MC-CRN xi init: all " << xi_start << std::endl;

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_MC_CRN);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_mc_crn(
            forwardGraph,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            xi_init,
            Arg._xi_lo,
            Arg._frim_mc,
            Arg._frim_max_sweeps);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_mc_crn";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_mc = Arg._frim_mc;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_naive,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM MC naive: benefit_inf_eval with fresh MC per J estimate, save xi.
    if (Arg._func == FRIM_MC_NAIVE)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        const double xi_start = Arg._xi_init_lo ? Arg._xi_lo : 1.0;
        std::vector<double> xi_init(graph.size(), xi_start);
        std::cout << "FRIM-MC-naive xi init: all " << xi_start << std::endl;

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_MC_NAIVE);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_mc_naive(
            forwardGraph,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            g_nodeHyperParams.q,
            xi_init,
            Arg._xi_lo,
            Arg._frim_mc,
            Arg._frim_max_sweeps,
            1e-9,
            Arg._casc_model);

        double q_sum = 0.0;
        for (double qv : g_nodeHyperParams.q)
            q_sum += qv;

        frimRes.run_info.method = "frim_mc_naive";
        frimRes.run_info.rand_seed = Arg._rand_seed;
        frimRes.run_info.num_v = graph.size();
        frimRes.run_info.num_mc = Arg._frim_mc;

        TIO::WriteFrimXi(Arg._outFileName, frimRes.xi, Arg._resultFolder);
        TIO::WriteFrimInfo(
            Arg._outFileName,
            frimRes.J_hat,
            frimRes.J_method_rr,
            frimRes.J_method_mc,
            frimRes.J_method_mc_naive,
            frimRes.J_method_rr_naive,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // Handle EVAL function: load xi, sample seeds from q, evaluate total benefit.
    if (Arg._func == EVAL)
    {
        const std::string evalResultFolder = Arg._resultFolder + "/evaluation";
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_for_eval;
        if (Arg._eval_xi_unirand_set || Arg._eval_xi_const_set)
        {
            if (Arg._eval_xi_unirand_set)
            {
                xi_for_eval.resize(graph.size());
                const double xi_hi = 1.0;
                const double span = xi_hi - Arg._xi_lo;
                for (size_t v = 0; v < graph.size(); v++)
                {
                    const double u = dsfmt_gv_genrand_close_open();
                    xi_for_eval[v] = Arg._xi_lo + u * span;
                }
                Arg.build_eval_xi_unirand_outfilename();
                std::cout << "Eval with xi ~ Uniform(" << Arg._xi_lo
                          << ", " << xi_hi << ") per node (unirand)" << std::endl;
            }
            else
            {
            xi_for_eval.assign(graph.size(), Arg._eval_xi_const);
            if (Arg._eval_xi_half_set)
            {
                const size_t half = graph.size() / 2;
                for (size_t v = 0; v < half; v++)
                    xi_for_eval[v] = Arg._xi_lo;
                for (size_t v = half; v < graph.size(); v++)
                    xi_for_eval[v] = 1.0;
                Arg.build_eval_xi_half_outfilename();
                std::cout << "Eval with xi: first " << half << " nodes = "
                          << Arg._xi_lo << ", rest = 1.0" << std::endl;
            }
            else
            {
                Arg.build_eval_xi_const_outfilename(Arg._eval_xi_const);
                std::cout << "Eval with xi ≡ " << Arg._eval_xi_const << " (constant)" << std::endl;
            }
            }
        }
        else
        {
            const FuncType frimMethod = Arg.decode_frim_eval_method();
            if (frimMethod != FRIM_RR && frimMethod != FRIM_RR_NAIVE
                && frimMethod != FRIM_RR_GRAPH
                && frimMethod != FRIM_RR_GRAPH_PRUNE
                && frimMethod != FRIM_SUBSIM
                && frimMethod != FRIM_MC_CRN && frimMethod != FRIM_MC_NAIVE)
            {
                std::cout << "Error: eval requires -method=frim_rr, -method=frim_rr_naive, "
                          << "-method=frim_rr_graph, "
                          << "-method=frim_rr_graph_prune, "
                          << "-method=frim_subsim, "
                          << "-method=frim_mc_crn, -method=frim_mc_naive, "
                          << "or -eval_xi=<value|half|unirand> for xi eval" << std::endl;
                return 0;
            }

            if (frimMethod == FRIM_RR_GRAPH || frimMethod == FRIM_RR_GRAPH_PRUNE)
            {
                const char* tag = (frimMethod == FRIM_RR_GRAPH)
                    ? "FRIM-RR-GRAPH eval"
                    : "FRIM-RR-GRAPH-PRUNE eval";
                double tau_sum = 0.0;
                for (double t : g_nodeHyperParams.tau) tau_sum += t;
                frimResolveRRGraphNumRR(Arg, graph.size(), tau_sum, tag);
            }

            Arg.build_frim_outfilename(frimMethod);
            if (!TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, g_frimXi, graph.size()))
            {
                std::cout << "Error: xi file not found for " << Arg._outFileName
                          << ", run the matching -func=frim_* solver first "
                          << "(use same -frim_rr or -eps, and -frim_sweeps)." << std::endl;
                return 0;
            }
            xi_for_eval = g_frimXi.xi;
        }

        const double res = GraphBase::benefit_inf_eval(
            graph,
            g_nodeHyperParams.q,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            Arg._casc_model,
            xi_for_eval,
            Arg._eval_mc);

        Graph rvsGraph;
        GraphBase::LoadGraph(rvsGraph, infilename);
        TResult rrRes;
        TAlg rrAlg(rvsGraph, rrRes);
        rrAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        const size_t num_rr = Arg._frim_rr;
        const size_t num_mc = Arg._eval_mc;
        // const double J_rr = rrAlg.frim_estimate_j_at_xi(
        //     xi_for_eval,
        //     g_nodeHyperParams.q,
        //     g_nodeHyperParams.tau,
        //     g_nodeHyperParams.lam,
        //     num_rr);
        // const double J_mc_gate = rrAlg.frim_estimate_mc_gate_at_xi(
        //     graph,
        //     xi_for_eval,
        //     g_nodeHyperParams.q,
        //     g_nodeHyperParams.tau,
        //     g_nodeHyperParams.lam,
        //     num_mc);

        if (!Arg._eval_xi_const_set && !Arg._eval_xi_unirand_set)
            std::cout << "Eval xi loaded from " << Arg._outFileName << std::endl;
        // std::cout << "  >>>RR J_hat (R=" << num_rr << "): " << J_rr << std::endl;
        // std::cout << "  >>>MC gate-objective (mc=" << num_mc << "): " << J_mc_gate << std::endl;
        if (Arg._eval_xi_const_set || Arg._eval_xi_unirand_set)
            std::cout << "  >>>MC benefit (eval, same mc): " << res << std::endl;
        std::cout << "  >>>MC-Benefit: " << res << std::endl;
        // std::cout << "  >>>RR vs MC-benefit gap: " << (J_rr - res)
        //           << " (" << (res > 0.0 ? 100.0 * (J_rr - res) / res : 0.0) << "%)" << std::endl;
        // std::cout << "  >>>RR vs MC-gate gap: " << (J_rr - J_mc_gate)
        //           << " (" << (J_mc_gate > 0.0 ? 100.0 * (J_rr - J_mc_gate) / J_mc_gate : 0.0)
        //           << "%)" << std::endl;

        // Fixed xi-independent baseline: sum of tau (full per-node value if not forwarding).
        // Do not use sum tau*(1-lam*xi): that moves with xi and breaks cross-xi ratio comparison.
        double total_node_benefit = 0.0;
        for (size_t v = 0; v < graph.size(); v++)
            total_node_benefit += g_nodeHyperParams.tau[v];

        TIO::WriteInfluence(Arg._outFileName, res, total_node_benefit, evalResultFolder);
        return 0;
    }

    // Initialize a result object to record the results
    TResult tRes;
    TAlg tAlg(graph, tRes);
    tAlg.set_vanilla_sample(Arg._vanilla);
    tAlg.set_prob_dist((ProbDist)probDist); // Set propagation model
    if (Arg._func == WIM) {
        tAlg.set_wrr_sample_mode(Arg._wrrSampleMode);
    }

    std::cout << "seedSize k=" << seedSize << std::endl;
    Arg.build_outfilename(seedSize, (ProbDist)probDist, graph, Arg._func);
    std::cout << "---The Begin of " << Arg._outFileName << "---\n";

    if (!Arg._hist)
    {
        if (Arg._func == IM) {
            if (Arg._num_samples > 0) {
                // Fixed number of RR sets sampling
                double actual_samples = Arg._num_samples;
                if (Arg._scaleN) {
                    actual_samples *= static_cast<double>(graph.size());
                }
                tAlg.fixed_subsim(seedSize, static_cast<int>(std::llround(actual_samples)));
            } else {
                tAlg.subsimOnly(seedSize, Arg._eps, delta);
            }
        } else if (Arg._func == WIM) {
            if (Arg._num_samples > 0) {
                double actual_samples = Arg._num_samples;
                if (Arg._scaleN) {
                    actual_samples *= static_cast<double>(graph.size());
                }
                tAlg.fixed_subsimW(seedSize, static_cast<int>(std::llround(actual_samples)));
            } else {
                tAlg.subsimWeight(seedSize, Arg._eps, delta);
            }
        }
    }
    else
    {
        std::cout <<"HIST is invoked." <<std::endl;
        if (seedSize < 10)
        {
            tAlg.subsimWithTrunc(seedSize, Arg._eps, delta);
        }
        else
        {
            tAlg.subsimWithHIST(seedSize, Arg._eps, delta);
        }
    }

    TIO::WriteResult(Arg._outFileName, tRes, Arg._resultFolder);
    TIO::WriteOrderSeeds(Arg._outFileName, tRes, Arg._resultFolder);
    std::cout << "---The End of " << Arg._outFileName << "---\n";
    tAlg.RefreshHypergraph();
    return 0;
}