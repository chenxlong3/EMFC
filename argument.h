#pragma once


class Argument
{
public:
    // Function parameter.
    // format: format graph
    // im: influence maximization
    std::string _funcStr = "im";
    FuncType _func = IM;

    // The number of nodes to be selected. Default is 50.
    int _seedsize = 50;

    // For the uniform setting, every edge has the same diffusion probability.
    float _probEdge = float(0.1);

    // Estimation error
    double _gamma = 0.05;

    // Error threshold 1-1/e-epsilon.
    double _eps = 0.1;

    // Failure probability delta. Default is 1/#nodes.
    double _delta = -1.0;

    // Graph name. Default is "facebook".
    std::string _graphname = "facebook";

    std::string _decay_mode = "exp";

    // Probability distribution
    // weights: graph data with weights
    // wc: wc setting
    // uniform: uniform setting
    // skewed: skewed distribution
    std::string _probDistStr = "wc";
    ProbDist _probDist = WC;
    CascadeModel _casc_model = IC;
    std::string _skewType = "exp";

    // Directory
    std::string _dir = "graphInfo";

    // Result folder
    std::string _resultFolder = "result";

    // File name of the result
    std::string _outFileName;

    std::string _candedges_filename, _seeds_filename;

    // wc variant
    double _wcVar = 1.0;

    // advanced or not
    bool _advanced = false;

    double _alpha = 0.5;
    
    //Random seed
    int _rand_seed = 0;
    // Number of samples for greedy
    double _num_samples = 0.0;
    // Scale num_samples by numV for actual RR sampling count
    bool _scaleN = false;
    // For evaluation
    std::string _method = "";

    // sample RR set with the vanilla method
    bool _vanilla = false;

    // Weighted RR root sampling for -func=wim: indeg|windeg|outdeg|uniform
    std::string _wrrSampleStr = "indeg";
    WRRSampleMode _wrrSampleMode = WRR_INDEG;

    // use hist algorithm
    bool _hist = false;

    // Node hyperparameter assignment (used by -func=format)
    std::string _lamModeStr = "two_tier";
    bool _lamModeError = false;
    double _lamUniformValue = 0.3;
    double _lamUnirandLo = 0.1;
    double _lamUnirandHi = 0.5;
    double _activeUserRatio = 0.2;
    double _activeQMean = 0.8;
    double _activeQVar = 1.0;
    double _inactiveQMean = 0.2;
    double _inactiveQVar = 1.0;

    Argument(int argc, char* argv[])
    {
        std::string param, value;

        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') break;

            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');

            if (!param.compare("-func")) _funcStr = value;
            else if (!param.compare("-seedsize")) _seedsize = stoi(value);
            else if (!param.compare("-num_samples")) _num_samples = stod(value);
            else if (!param.compare("-eps")) _eps = stod(value);
            else if (!param.compare("-delta")) _delta = stod(value);
            else if (!param.compare("-gamma")) _gamma = stod(value);
            else if (!param.compare("-model")) _casc_model = value == "LT" ? LT : IC;
            else if (!param.compare("-gname")) _graphname = value;
            else if (!param.compare("-dir")) _dir = value;
            else if (!param.compare("-outpath")) _resultFolder = value;
            else if (!param.compare("-pdist")) _probDistStr = value;
            else if (!param.compare("-pedge")) _probEdge = stof(value);
            else if (!param.compare("-wcvariant")) _wcVar = stod(value);
            else if (!param.compare("-skew")) _skewType = value;
            else if (!param.compare("-vanilla")) _vanilla = (value == "1");
            else if (!param.compare("-hist")) _hist = (value == "1");
            else if (!param.compare("-advanced")) _advanced = (value == "1");
            else if (!param.compare("-rand_seed")) _rand_seed = stoi(value);
            else if (!param.compare("-method")) _method = value;
            else if (!param.compare("-alpha")) _alpha = stof(value);
            else if (!param.compare("-decay_mode")) _decay_mode = value;
            else if (!param.compare("-scaleN")) _scaleN = (value == "1");
            else if (!param.compare("-wrr_sample")) _wrrSampleStr = value;
            else if (!param.compare("-lam_mode")) _lamModeStr = value;
            else if (!param.compare("-lam_value")) _lamUniformValue = stod(value);
            else if (!param.compare("-lam_lo")) _lamUnirandLo = stod(value);
            else if (!param.compare("-lam_hi")) _lamUnirandHi = stod(value);
            else if (!param.compare("-active_ratio")) _activeUserRatio = stod(value);
            else if (!param.compare("-active_q_mean")) _activeQMean = stod(value);
            else if (!param.compare("-active_q_var")) _activeQVar = stod(value);
            else if (!param.compare("-inactive_q_mean")) _inactiveQMean = stod(value);
            else if (!param.compare("-inactive_q_var")) _inactiveQVar = stod(value);
        }

        if (_wcVar <= 0)
        {
            //wrong input
            _wcVar = 1.0;
        }

        decode_wrr_sample_mode();
        decode_lam_mode();

        // Fixed-sample mode: epsilon should be disabled.
        if (_num_samples > 0.0)
        {
            _eps = 0.0;
        }

        decode_func_type();
        decode_prob_dist();
    }



    void build_outfilename(int seedSize, ProbDist dist, Graph& graph, FuncType func=IM, std::string method="ADIL")
    {
        std::string distStr; 

        if (dist == WEIGHTS)
        {
            _probDistStr = "weights";
        }
        else if (dist == WC)
        {
            _probDistStr = "wc";
        }
        else if (dist == UNIFORM)
        {
            _probDistStr = "uniform";

            for (int i = 0; i < graph.size(); i++)
            {
                if (graph[i].size() > 0 )
                {
                    _probEdge = graph[i][0].second;
                    break;
                }
            }
        }
        else
        {
            _probDistStr = "skewed";
        }
        if (_func == IM) {
            std::string method = "subsim";
            if (this->_num_samples > 0.0) method = method + "[" + format_num_samples() + "]";
            _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
        }
        else if (_func == WIM) {
            std::string method = wrr_method_tag();
            if (this->_num_samples > 0.0) method = method + "[" + format_num_samples() + "]";
            _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
        }
        else if (func == ScaLIM) {
            std::string method = "ScaLIM";
            if (this->_num_samples >= 1) method = method + "[" + format_num_samples() + "]";
            _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
            _outFileName += "_" + this->_decay_mode + std::to_string(this->_alpha);
        }
        else if (func == RAND) {
            _outFileName = TIO::BuildOutFileName(_graphname, "RAND", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0);
        }
        else if (func == OUTDEG) {
            _outFileName = TIO::BuildOutFileName(_graphname, "OUTDEG", seedSize, _probDistStr, _probEdge, 0, 0.0);
        }
        else if (func == PROB) {
            _outFileName = TIO::BuildOutFileName(_graphname, "PROB", seedSize, _probDistStr, _probEdge, 0, 0.0);
        }
        else if (func == DG) {
            _outFileName = TIO::BuildOutFileName(_graphname, "DG", seedSize, _probDistStr, _probEdge, 0, 0.0);
            _outFileName += "_" + this->_decay_mode + std::to_string(this->_alpha);
        }
        else if (func == GREEDY) {
            _outFileName = TIO::BuildOutFileName(_graphname, "GREEDY[" + format_num_samples() + "]", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0);
            _outFileName += "_" + this->_decay_mode + std::to_string(this->_alpha);
        }
        else if (func == SUBSIM)
        {
            std::string method = "SUBSIM";
            if (this->_num_samples >= 1) method = method + "[" + format_num_samples() + "]";
            _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
            _outFileName += "_" + this->_decay_mode + std::to_string(this->_alpha);
        }
        else if (func == AIS)
        {
            _outFileName = TIO::BuildOutFileName(_graphname, "AIS", seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
        }
        else if (func == EVAL)
        {
            if(this->_method == "subsim") {
                std::string method = this->_method;
                if (this->_num_samples > 0.0) method = method + "[" + format_num_samples() + "]";
                _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
            }
            else if (this->_method == "subsimW" || this->_method == "subsimW-win"
                     || this->_method == "subsimW-wout" || this->_method == "subsimW-out"
                     || this->_method == "subsimW-uni") {
                std::string method = this->_method;
                if (this->_num_samples > 0.0) method = method + "[" + format_num_samples() + "]";
                _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
            }
            else if (this->_method == "GREEDY") {
                _outFileName = TIO::BuildOutFileName(_graphname, "GREEDY[" + format_num_samples() + "]", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0);
            }
            else if (this->_method == "OUTDEG" || this->_method == "PROB" || this->_method == "DG")
            {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, 0, 0.0);
            }
            else if (this->_method == "AIS") {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0);
            }
            else if (this->_method == "ADIL" || this->_method == "ScaLIM" || this->_method == "SUBSIM")
            {
                std::string method = this->_method;
                if (this->_num_samples >= 0.0) method = method + "[" + format_num_samples() + "]";
                _outFileName = TIO::BuildOutFileName(_graphname, method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
            }
            
        }

        return ;
    }

    // Fill candedges and seeds filenames from Arg
    void build_cand_seeds_filenames() {
        this->_candedges_filename = this->_dir + "/" + "candEdges_" + this->_graphname + "_num" + std::to_string(this->_seedsize);
        this->_seeds_filename = this->_dir + "/seed/" + "seed_" + this->_graphname + "_num" + std::to_string(this->_seedsize);
    }

    void decode_prob_dist()
    {
        if (_probDistStr == "wc")
        {
            _probDist = WC;
        }
        else if (_probDistStr == "wc2")
        {
            _probDist = WC2;
        }
        else if (_probDistStr == "uniform")
        {
            _probDist = UNIFORM;
        }
        else if (_probDistStr == "skewed")
        {
            _probDist = SKEWED;
        }
        else if (_probDistStr == "weights")
        {
            _probDist = WEIGHTS;
        }
        else 
        {
            _probDist = PROB_DIST_ERROR;
        }
    }

    void decode_func_type()
    {
        if (_funcStr == "format")
        {
            _func = FORMAT;
        }
        else if (_funcStr == "im")
        {
            _func = IM;
        }
        else if (_funcStr == "wim")
        {
            _func = WIM;
        }
        else if (_funcStr == "SUBSIM")
        {
            _func = SUBSIM;
        }
        else if (_funcStr == "OUTDEG")
        {
            _func = OUTDEG;
        }
        // PROB
        else if (_funcStr == "PROB")
        {
            _func = PROB;
        }
        // RAND
        else if (_funcStr == "RAND")
        {
            _func = RAND;
        }
        // AIS
        else if (_funcStr == "AIS")
        {
            _func = AIS;
        }
        // SANDWICH
        else if (_funcStr == "ScaLIM")
        {
            _func = ScaLIM;
        }
        // ADIL
        else if (_funcStr == "ADIL")
        {
            _func = ADIL;
        }
        // GREEDY
        else if (_funcStr == "GREEDY")
        {
            _func = GREEDY;
        }
        else if (_funcStr == "CHANGE_PROB")
        {
            _func = CHANGE_PROB;
        }
        else if (_funcStr == "DG")
        {
            _func = DG;
        }
        else if (_funcStr == "EVAL" || _funcStr == "eval")
        {
            _func = EVAL;
        }
        else if (_funcStr == "PREP_CAND")
        {
            _func = PREP_CAND;
        }
        else if (_funcStr == "STAT" || _funcStr == "stat")
        {
            _func = STAT;
        }
        else
        {
            _func = FUNC_ERROR;
        }
    }

    std::string format_num_samples() const
    {
        if (std::floor(_num_samples) == _num_samples)
        {
            return std::to_string(static_cast<long long>(_num_samples));
        }
        std::ostringstream oss;
        oss << _num_samples;
        return oss.str();
    }

    std::string wrr_method_tag() const
    {
        switch (_wrrSampleMode) {
        case WRR_WINDEG: return "subsimW-win";
        case WRR_WOUTDEG: return "subsimW-wout";
        case WRR_OUTDEG: return "subsimW-out";
        case WRR_UNIFORM: return "subsimW-uni";
        case WRR_INDEG:
        default: return "subsimW";
        }
    }

    void decode_lam_mode()
    {
        if (_lamModeStr == "uniform" || _lamModeStr == "unirand" || _lamModeStr == "two_tier")
        {
            _lamModeError = false;
        }
        else
        {
            _lamModeError = true;
        }
    }

    void decode_wrr_sample_mode()
    {
        if (_wrrSampleStr == "windeg" || _wrrSampleStr == "w_in"
            || _wrrSampleStr == "weighted_indegree" || _wrrSampleStr == "weighted-indegree") {
            _wrrSampleMode = WRR_WINDEG;
        } else if (_wrrSampleStr == "woutdeg" || _wrrSampleStr == "w_out"
                   || _wrrSampleStr == "weighted_outdegree" || _wrrSampleStr == "weighted-outdegree") {
            _wrrSampleMode = WRR_WOUTDEG;
        } else if (_wrrSampleStr == "outdeg" || _wrrSampleStr == "out"
                   || _wrrSampleStr == "outdegree" || _wrrSampleStr == "out-degree") {
            _wrrSampleMode = WRR_OUTDEG;
        } else if (_wrrSampleStr == "uniform" || _wrrSampleStr == "uni") {
            _wrrSampleMode = WRR_UNIFORM;
        } else if (_wrrSampleStr == "indeg" || _wrrSampleStr == "in"
                   || _wrrSampleStr == "indegree" || _wrrSampleStr == "in-degree") {
            _wrrSampleMode = WRR_INDEG;
        } else {
            _wrrSampleMode = WRR_SAMPLE_ERROR;
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
