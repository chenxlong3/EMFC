#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

class IoController
{
public:
    static void mkdir_absence(const char* outFolder)
    {
#if defined(_WIN32)
        CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
        mkdir(outFolder, 0733); // can be used on non-Windows
#endif
    }

    /// Save a serialized file
    template <class T>
    static void SaveFile(const std::string filename, const T& output)
    {
        std::ofstream outfile(filename, std::ios::binary);

        if (!outfile.eof() && !outfile.fail())
        {
            StreamType res;
            serialize(output, res);
            outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
            outfile.close();
            res.clear();
            std::cout << "Save file successfully: " << filename << '\n';
        }
        else
        {
            std::cout << "Save file failed: " + filename << '\n';
            exit(1);
        }
    }

    /// Load a serialized file
    template <class T>
    static void LoadFile(const std::string filename, T& input)
    {
        std::ifstream infile(filename, std::ios::binary);

        if (!infile.eof() && !infile.fail())
        {
            infile.seekg(0, std::ios_base::end);
            const std::streampos fileSize = infile.tellg();
            infile.seekg(0, std::ios_base::beg);
            std::vector<uint8_t> res(fileSize);
            infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
            infile.close();
            input.clear();
            auto it = res.cbegin();
            input = deserialize<T>(it, res.cend());
            res.clear();
        }
        else
        {
            std::cout << "Cannot open file: " + filename << '\n';
            exit(1);
        }
    }

    /// Save graph structure to a file
    static void SaveGraphStruct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        SaveFile(filename, vecGraph);
    }

    /// Load graph structure from a file
    static void LoadGraphStruct(const std::string graphName, Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        LoadFile(filename, vecGraph);
    }

    static void SaveGraphProbDist(const std::string graphName, int dist)
    {
        std::ofstream outFile(graphName + ".probdist");
        outFile<< dist;
    }


    static int LoadGraphProbDist(const std::string graphName)
    {
        std::string filename = graphName + ".probdist";
        std::ifstream infile(filename);
        int probDist = WEIGHTS;

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return probDist;
        }

        infile >> probDist;
        infile.close();
        std::cout << "probability distribution: " << probDist << std::endl;
        return probDist;
    }

    static void SaveGraphNodeHyperParams(
        const std::string& graphName,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam)
    {
        const std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> blob(q, tau, lam);
        SaveFile(graphName + ".nodehyper.vec", blob);
    }

    static bool readGraphNodeHyperParamsFile(
        const std::string& graphName,
        std::vector<double>& q,
        std::vector<double>& tau,
        std::vector<double>& lam)
    {
        const std::string filename = graphName + ".nodehyper.vec";
        std::ifstream infile(filename, std::ios::binary);
        if (infile.eof() || infile.fail())
        {
            return false;
        }

        infile.seekg(0, std::ios_base::end);
        const std::streampos fileSize = infile.tellg();
        infile.seekg(0, std::ios_base::beg);
        std::vector<uint8_t> res(static_cast<size_t>(fileSize));
        infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
        infile.close();

        auto it = res.cbegin();
        std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> blob =
            deserialize<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>(it, res.cend());
        q = std::move(std::get<0>(blob));
        tau = std::move(std::get<1>(blob));
        lam = std::move(std::get<2>(blob));
        return true;
    }

    static void LoadGraphNodeHyperParams(
        const std::string& graphName,
        std::vector<double>& q,
        std::vector<double>& tau,
        std::vector<double>& lam)
    {
        if (!readGraphNodeHyperParamsFile(graphName, q, tau, lam))
        {
            std::cout << "Cannot open file: " + graphName + ".nodehyper.vec" << '\n';
            exit(1);
        }
    }

    /// Load q/tau/lam into global-style storage; returns false if file is missing.
    static bool LoadGraphNodeHyperParamsData(const std::string& graphName, NodeHyperParamsData& out)
    {
        out.q.clear();
        out.tau.clear();
        out.lam.clear();
        out.loaded = false;

        if (!readGraphNodeHyperParamsFile(graphName, out.q, out.tau, out.lam))
        {
            return false;
        }

        out.loaded = true;
        std::cout << "Loaded node hyper params: q/tau/lam size = "
                  << out.q.size() << std::endl;
        return true;
    }

    /// Get out-file name
    static std::string BuildOutFileName(const std::string graphName, const std::string algName, const int seedsize,
                                        const std::string probDist, const float probEdge, const int rand_seed = 0, const double eps = 0.1)
    {
        if (probDist == "uniform")
        {
            return graphName + "_" + std::to_string(rand_seed) + "_" + algName + "_k" + std::to_string(seedsize) + "_" + std::to_string(eps)  + "_" + probDist + "_" + std::
                   to_string(probEdge);
        }

        return graphName + "_" + std::to_string(rand_seed) + "_" + algName + "_k" + std::to_string(seedsize) + "_" + std::to_string(eps) + "_" + probDist;
    }

    /// Print the results
    static void WriteResult(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        const auto approx = resultObj.get_approximation();
        const auto runTime = resultObj.get_running_time();
        const auto influence = resultObj.get_influence();
        const auto influenceOriginal = resultObj.get_influence_original();
        const auto seedSize = resultObj.get_seed_size();
        const auto RRsetsSize = resultObj.get_RRsets_size();
        const auto samplingTime = resultObj.get_sampling_time();
        const auto selectionTime = resultObj.get_selection_time();
        std::cout << "   --------------------" << std::endl;
        std::cout << "  |Approx.: " << approx << std::endl;
        std::cout << "  |Time (sec): " << runTime << std::endl;
        std::cout << "  |Sampling Time (sec): " << samplingTime << std::endl;
        std::cout << "  |Selection Time (sec): " << selectionTime << std::endl;
        std::cout << "  |Influence: " << influence << std::endl;
        std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
        std::cout << "  |#Seeds: " << seedSize << std::endl;
        std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
        std::cout << "   --------------------" << std::endl;
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/info";
        mkdir_absence(outpath.c_str());
        std::ofstream outFileNew(outpath + "/" + outFileName);

        if (outFileNew.is_open())
        {
            outFileNew << "Approx.: " << approx << std::endl;
            outFileNew << "Time (sec): " << runTime << std::endl;
            outFileNew << "Sampling Time (sec): " << samplingTime << std::endl;
            outFileNew << "Selection Time (sec): " << selectionTime << std::endl;
            outFileNew << "Influence: " << influence << std::endl;
            outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
            outFileNew << "#Seeds: " << seedSize << std::endl;
            outFileNew << "#RR sets: " << RRsetsSize << std::endl;
            outFileNew.close();
        }
    }

    /// Print the seeds
    static void WriteOrderSeeds(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        auto vecSeed = resultObj.get_seed_vec();
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/seed";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile(outpath + "/seed_" + outFileName);

        for (auto i = 0; i < vecSeed.size(); i++)
        {
            outFile << vecSeed[i] << '\n';
        }

        outFile.close();
    }

    /// Read edges from result file (for evaluation)
    static void ReadEdges(const std::string& outFileName, std::vector<std::pair<uint32_t, uint32_t>>& edges, const std::string& resultFolder)
    {
        edges.clear();
        std::ifstream inFile(resultFolder + "/edge/" + outFileName);
        if (!inFile.is_open())
        {
            std::cout << "Cannot open edge file: " << resultFolder + "/edge/" + outFileName << std::endl;
            return;
        }
        uint32_t u, v;
        while (inFile >> u >> v)
        {
            edges.push_back(std::make_pair(u, v));
        }
        inFile.close();
    }

    /// Read seeds from seed file (for evaluation)
    static void ReadSeeds(const std::string& seedsFile, std::vector<uint32_t>& seeds)
    {
        seeds.clear();
        std::ifstream inFile(seedsFile);
        if (!inFile.is_open())
        {
            std::cout << "Cannot open seeds file: " << seedsFile << std::endl;
            return;
        }
        uint32_t seed;
        while (inFile >> seed)
        {
            seeds.push_back(seed);
        }
        inFile.close();
    }

    /// Write influence result (for evaluation)
    static void WriteInfluence(
        const std::string& outFileName,
        const double influence,
        const double totalNodeBenefit,
        const std::string& resultFolder)
    {
        mkdir_absence(resultFolder.c_str());
        std::ofstream outFile(resultFolder + "/eval_" + outFileName);

        if (outFile.is_open())
        {
            const double ratio = (totalNodeBenefit <= 0.0)
                ? 0.0
                : (influence / totalNodeBenefit);
            outFile << "Influence: " << influence << std::endl;
            outFile << "TotalNodeBenefit: " << totalNodeBenefit << std::endl;
            outFile << "Ratio: " << ratio << std::endl;
            outFile.close();
            std::cout << "Evaluation result saved to: " << resultFolder + "/eval_" + outFileName << std::endl;
        }
        else
        {
            std::cout << "Failed to write evaluation result" << std::endl;
        }
    }

    /// Save FRIM xi vector (one value per line), mirroring WriteOrderSeeds layout.
    static void WriteFrimXi(
        const std::string& outFileName,
        const std::vector<double>& xi,
        const std::string& resultFolder)
    {
        mkdir_absence(resultFolder.c_str());
        const auto outpath = resultFolder + "/xi";
        mkdir_absence(outpath.c_str());
        const std::string filepath = outpath + "/xi_" + outFileName;
        std::ofstream outFile(filepath);

        for (size_t i = 0; i < xi.size(); i++)
        {
            outFile << xi[i] << '\n';
        }

        outFile.close();
        std::cout << "Saved xi to: " << filepath << std::endl;
    }

    static bool readFrimXiFile(
        const std::string& resultFolder,
        const std::string& outFileName,
        std::vector<double>& xi)
    {
        const std::string filepath = resultFolder + "/xi/xi_" + outFileName;
        std::ifstream inFile(filepath);
        if (!inFile.is_open())
        {
            return false;
        }

        xi.clear();
        double val = 0.0;
        while (inFile >> val)
        {
            xi.push_back(val);
        }
        inFile.close();
        return true;
    }

    /// Load xi into global-style storage; returns false if file is missing or size mismatches.
    static bool LoadFrimXiData(
        const std::string& resultFolder,
        const std::string& outFileName,
        FrimXiData& out,
        size_t expectedSize)
    {
        out.xi.clear();
        out.J_hat = 0.0;
        out.J_method_rr = 0.0;
        out.J_method_mc = 0.0;
        out.loaded = false;

        if (!readFrimXiFile(resultFolder, outFileName, out.xi))
        {
            return false;
        }

        if (out.xi.size() != expectedSize)
        {
            std::cout << "Error: xi size " << out.xi.size()
                      << " != graph size " << expectedSize << std::endl;
            return false;
        }

        const std::string infoPath = resultFolder + "/info/" + outFileName;
        std::ifstream infoFile(infoPath);
        if (infoFile.is_open())
        {
            std::string line;
            while (std::getline(infoFile, line))
            {
                if (line.find("J_hat:") == 0)
                    out.J_hat = std::stod(line.substr(7));
                else if (line.find("J_method_rr:") == 0)
                    out.J_method_rr = std::stod(line.substr(13));
                else if (line.find("J_method_mc:") == 0)
                    out.J_method_mc = std::stod(line.substr(13));
            }
            infoFile.close();
        }

        out.loaded = true;
        std::cout << "Loaded xi: size = " << out.xi.size() << std::endl;
        return true;
    }

    /// Save FRIM summary info, mirroring WriteResult layout.
    static void WriteFrimInfo(
        const std::string& outFileName,
        double J_hat,
        double J_method_rr,
        double J_method_mc,
        double J_method_mc_naive,
        double J_method_rr_naive,
        size_t seedCount,
        const FrimRunInfo& runInfo,
        const std::string& resultFolder)
    {
        std::cout << "   --------------------" << std::endl;
        std::cout << "  |J_hat: " << J_hat << std::endl;
        std::cout << "  |J_method_rr: " << J_method_rr << std::endl;
        std::cout << "  |J_method_rr_naive: " << J_method_rr_naive << std::endl;
        std::cout << "  |J_method_mc: " << J_method_mc << std::endl;
        std::cout << "  |J_method_mc_naive: " << J_method_mc_naive << std::endl;
        std::cout << "  |E[|seeds|] sum q: " << seedCount << std::endl;
        std::cout << "  |method: " << runInfo.method << std::endl;
        std::cout << "  |time_total_sec: " << runInfo.time_total_sec << std::endl;
        std::cout << "  |time_sample_sec: " << runInfo.time_sample_sec << std::endl;
        std::cout << "  |time_solve_sec: " << runInfo.time_solve_sec << std::endl;
        std::cout << "  |time_estimate_sec: " << runInfo.time_estimate_sec << std::endl;
        std::cout << "   --------------------" << std::endl;

        mkdir_absence(resultFolder.c_str());
        const auto outpath = resultFolder + "/info";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile(outpath + "/" + outFileName);

        if (outFile.is_open())
        {
            outFile << "J_hat: " << J_hat << std::endl;
            outFile << "J_method_rr: " << J_method_rr << std::endl;
            outFile << "J_method_rr_naive: " << J_method_rr_naive << std::endl;
            outFile << "J_method_mc: " << J_method_mc << std::endl;
            outFile << "J_method_mc_naive: " << J_method_mc_naive << std::endl;
            outFile << "E[|seeds|] sum q: " << seedCount << std::endl;
            outFile << "method: " << runInfo.method << std::endl;
            outFile << "rand_seed: " << runInfo.rand_seed << std::endl;
            outFile << "num_v: " << runInfo.num_v << std::endl;
            outFile << "num_rr: " << runInfo.num_rr << std::endl;
            outFile << "num_mc: " << runInfo.num_mc << std::endl;
            outFile << "max_sweeps: " << runInfo.max_sweeps << std::endl;
            outFile << "sweeps_completed: " << runInfo.sweeps_completed << std::endl;
            outFile << "xi_lo: " << runInfo.xi_lo << std::endl;
            outFile << "num_xi_one: " << runInfo.num_xi_one << std::endl;
            outFile << "num_xi_lo: " << runInfo.num_xi_lo << std::endl;
            outFile << "time_total_sec: " << runInfo.time_total_sec << std::endl;
            outFile << "time_sample_sec: " << runInfo.time_sample_sec << std::endl;
            outFile << "time_solve_sec: " << runInfo.time_solve_sec << std::endl;
            outFile << "time_estimate_sec: " << runInfo.time_estimate_sec << std::endl;
            outFile.close();
            std::cout << "Saved FRIM info to: " << outpath + "/" + outFileName << std::endl;
        }
    }
};

using TIO = IoController;
using PIO = std::shared_ptr<IoController>;
