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
    static void WriteInfluence(const std::string& outFileName, const double influence, const size_t numV, const std::string& resultFolder)
    {
        mkdir_absence(resultFolder.c_str());
        std::ofstream outFile(resultFolder + "/eval_" + outFileName);

        if (outFile.is_open())
        {
            const double ratio = (numV == 0) ? 0.0 : (influence / static_cast<double>(numV));
            outFile << "Influence: " << influence << std::endl;
            outFile << "Ratio: " << ratio << std::endl;
            outFile.close();
            std::cout << "Evaluation result saved to: " << resultFolder + "/eval_" + outFileName << std::endl;
        }
        else
        {
            std::cout << "Failed to write evaluation result" << std::endl;
        }
    }
};

using TIO = IoController;
using PIO = std::shared_ptr<IoController>;
