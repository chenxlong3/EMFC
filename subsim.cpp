#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

void init_random_seed()
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
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
        LogInfo("The supported func: format, im, wim, SUBSIM, OUTDEG, PROB, RAND, AIS, ScaLIM, ADIL, GREEDY, DG, EVAL, PREP_CAND, STAT");
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

    init_random_seed();

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
    // Load graph: evaluation needs forward graph, optimization uses reverse graph.
    Graph graph;
    if (Arg._func == EVAL) {
        GraphBase::LoadForwardGraph(graph, infilename);
    } else {
        GraphBase::LoadGraph(graph, infilename);
    }
    int probDist = GraphBase::LoadGraphProbDist(infilename);

    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();

    // Handle EVAL function: sample seeds from q, evaluate total benefit.
    if (Arg._func == EVAL)
    {
        const std::string evalResultFolder = "./result/evaluation";
        if (!g_nodeHyperParams.loaded)
        {
            std::cout << "Error: node hyper params not loaded for " << infilename
                      << ", run -func=format first." << std::endl;
            return 0;
        }

        std::vector<double> xi(graph.size(), 0.9);

        const double res = GraphBase::benefit_inf_eval(
            graph,
            g_nodeHyperParams.q,
            g_nodeHyperParams.tau,
            g_nodeHyperParams.lam,
            Arg._casc_model,
            xi);

        const std::string outFileName = Arg._graphname + "_" + Arg._probDistStr + "_benefit";
        TIO::WriteInfluence(outFileName, res, graph.size(), evalResultFolder);
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