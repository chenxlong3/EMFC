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
        LogInfo("The supported func: format, im, wim, SUBSIM, OUTDEG, PROB, RAND, AIS, ScaLIM, ADIL, GREEDY, DG, EVAL, PREP_CAND, STAT, frim_rr, frim_rr_naive, frim_rr_crn, frim_mc, frim_mc_naive");
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
        LogInfo("The supported lam_mode: uniform, unirand, two_tier");
        return 0;
    }

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
                Arg._inactiveQVar));
        return 0;
    }

    if (!TIO::LoadGraphNodeHyperParamsData(infilename, g_nodeHyperParams))
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
    else if (Arg._func == FRIM_MC || Arg._func == FRIM_MC_NAIVE)
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

        std::vector<double> xi_init(graph.size(), 1.0);
        Arg.build_frim_outfilename(FRIM_RR);
        FrimXiData xiRR;
        if (TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, xiRR, graph.size()))
        {
            xi_init = xiRR.xi;
            std::cout << "FRIM-RR-naive xi init loaded from frim_rr: " << Arg._outFileName << std::endl;
        }
        else
        {
            std::cout << "FRIM-RR-naive xi init: all 1.0 (frim_rr result not found)" << std::endl;
        }

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

    // FRIM RR CRN: fixed tau-root structural samples, xi_gate CRN per node.
    if (Arg._func == FRIM_RR_CRN)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_init(graph.size(), 1.0);
        Arg.build_frim_outfilename(FRIM_RR);
        FrimXiData xiRR;
        if (TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, xiRR, graph.size()))
        {
            xi_init = xiRR.xi;
            std::cout << "FRIM-RR-CRN xi init loaded from frim_rr: " << Arg._outFileName << std::endl;
        }
        else
        {
            std::cout << "FRIM-RR-CRN xi init: all 1.0 (frim_rr result not found)" << std::endl;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_RR_CRN);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_rr_crn(
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

        frimRes.run_info.method = "frim_rr_crn";
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
            frimRes.J_method_rr_crn,
            static_cast<size_t>(std::llround(q_sum)),
            frimRes.run_info,
            Arg._resultFolder);
        std::cout << "---The End of " << Arg._outFileName << "---\n";
        return 0;
    }

    // FRIM MC method: CRN coordinate ascent, save xi.
    if (Arg._func == FRIM_MC)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_init(graph.size(), 1.0);
        Arg.build_frim_outfilename(FRIM_RR);
        FrimXiData xiRR;
        if (TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, xiRR, graph.size()))
        {
            xi_init = xiRR.xi;
            std::cout << "FRIM-MC xi init loaded from frim_rr: " << Arg._outFileName << std::endl;
        }
        else
        {
            std::cout << "FRIM-MC xi init: all 1.0 (frim_rr result not found)" << std::endl;
        }

        TResult tRes;
        TAlg tAlg(graph, tRes);
        tAlg.set_prob_dist(static_cast<ProbDist>(probDist));

        Arg.build_frim_outfilename(FRIM_MC);
        std::cout << "---The Begin of " << Arg._outFileName << "---\n";

        FrimXiResult frimRes = tAlg.frim_solve_mc(
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

        frimRes.run_info.method = "frim_mc";
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

    // FRIM MC naive: full BFS delta on CRN live-edge samples, save xi.
    if (Arg._func == FRIM_MC_NAIVE)
    {
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_init(graph.size(), 1.0);
        Arg.build_frim_outfilename(FRIM_RR);
        FrimXiData xiRR;
        if (TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, xiRR, graph.size()))
        {
            xi_init = xiRR.xi;
            std::cout << "FRIM-MC-naive xi init loaded from frim_rr: " << Arg._outFileName << std::endl;
        }
        else
        {
            std::cout << "FRIM-MC-naive xi init: all 1.0 (frim_rr result not found)" << std::endl;
        }

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
            Arg._frim_max_sweeps);

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
        const std::string evalResultFolder = "./result/evaluation";
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi_for_eval;
        if (Arg._eval_xi_const_set)
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
        else
        {
            const FuncType frimMethod = Arg.decode_frim_eval_method();
            if (frimMethod != FRIM_RR && frimMethod != FRIM_RR_NAIVE
                && frimMethod != FRIM_RR_CRN
                && frimMethod != FRIM_MC && frimMethod != FRIM_MC_NAIVE)
            {
                std::cout << "Error: eval requires -method=frim_rr, -method=frim_rr_naive, "
                          << "-method=frim_rr_crn, "
                          << "-method=frim_mc, -method=frim_mc_naive, "
                          << "or -eval_xi=<value> for constant xi" << std::endl;
                return 0;
            }

            Arg.build_frim_outfilename(frimMethod);
            if (!TIO::LoadFrimXiData(Arg._resultFolder, Arg._outFileName, g_frimXi, graph.size()))
            {
                std::cout << "Error: xi file not found for " << Arg._outFileName
                          << ", run the matching -func=frim_* solver first." << std::endl;
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
        const double J_rr = rrAlg.frim_estimate_j_at_xi(
            xi_for_eval,
            g_nodeHyperParams.q,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            num_rr);
        const double J_mc_gate = rrAlg.frim_estimate_mc_gate_at_xi(
            graph,
            xi_for_eval,
            g_nodeHyperParams.q,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            num_mc);

        if (!Arg._eval_xi_const_set)
            std::cout << "Eval xi loaded from " << Arg._outFileName << std::endl;
        std::cout << "  >>>RR J_hat (R=" << num_rr << "): " << J_rr << std::endl;
        std::cout << "  >>>MC gate-objective (mc=" << num_mc << "): " << J_mc_gate << std::endl;
        if (Arg._eval_xi_const_set)
            std::cout << "  >>>MC benefit (eval, same mc): " << res << std::endl;
        std::cout << "  >>>RR vs MC-benefit gap: " << (J_rr - res)
                  << " (" << (res > 0.0 ? 100.0 * (J_rr - res) / res : 0.0) << "%)" << std::endl;
        std::cout << "  >>>RR vs MC-gate gap: " << (J_rr - J_mc_gate)
                  << " (" << (J_mc_gate > 0.0 ? 100.0 * (J_rr - J_mc_gate) / J_mc_gate : 0.0)
                  << "%)" << std::endl;

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