#pragma once

#include <cmath>
#include <algorithm>

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
    double _eps = 0.5;
    /// Set when -eps is passed on the command line.
    bool _eps_set = false;

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
    bool _eval_xi_const_set = false;
    bool _eval_xi_half_set = false;
    bool _eval_xi_unirand_set = false;
    double _eval_xi_const = 1.0;
    uint32_t _eval_mc = 10000;

    // sample RR set with the vanilla method
    bool _vanilla = false;

    // Weighted RR root sampling for -func=wim: indeg|windeg|outdeg|uniform
    std::string _wrrSampleStr = "indeg";
    WRRSampleMode _wrrSampleMode = WRR_INDEG;

    // use hist algorithm
    bool _hist = false;

    // Node hyperparameter assignment (used by -func=format)
    std::string _lamModeStr = "unirand";
    bool _lamModeError = false;
    double _lamUniformValue = 0.3;
    double _lamUnirandLo = 0.1;
    double _lamUnirandHi = 0.5;
    double _activeUserRatio = 0.2;
    double _activeQMean = 0.8;
    double _activeQVar = 1.0;
    double _inactiveQMean = 0.2;
    double _inactiveQVar = 1.0;

    std::string _qModeStr = "exponential";
    bool _qModeError = false;
    double _qExpScale = 1.0;
    double _qUnirandLo = 0.0;
    double _qUnirandHi = 1.0;

    /// Optional suffix for loading {graph}.nodehyper.{suffix}.vec at runtime.
    std::string _nodehyperSuffix;
    /// Base suffix when deriving from an existing nodehyper file (empty = default).
    std::string _nodehyperBaseSuffix;
    /// Output suffix for -func=hyp_derive (also used as load suffix if _nodehyperSuffix empty).
    std::string _hypOutputSuffix;
    std::string _hypProfilePath;

    std::string _tauModeStr = "q_normal";
    bool _tauModeError = false;
    double _tauLo = 1.0;
    double _tauHi = 5.0;
    double _tauJitter = -1.0;
    double _tauQVar = 1.0;

    // FRIM xi selection
    double _xi_lo = 0.5;
    // rr-graph constant init: false=all 1.0, true=all xi_lo (-xi_init=lo)
    bool _xi_init_lo = false;
    size_t _frim_rr = 10000;
    /// Set when -frim_rr is passed; otherwise rr-graph uses the delta-based formula in resolveFrimRRGraphNumSamples.
    bool _frim_rr_fixed = false;
    size_t _frim_mc = 1000;
    int _frim_max_sweeps = 10;
    /// rr-graph sweep index: omit (u,r) when xi[u] in {xi_lo,1} cannot change sample r's weight.
    bool _rr_graph_gate_sweep_index = true;
    /// When true (default): discard no-hit RR-graphs after sampling; R still counts all attempts.
    /// Pass -rr_graph_keep_no_hit=1 to store all samples (legacy behavior).
    bool _rr_graph_store_hit_only = true;
    /// frim_rr_root_stat: also build RR-graphs (slow) to count stored roots.
    bool _rr_root_stat_build = false;
    /// rr-graph warm-start: sample first R/div graphs, run warm sweeps, then sample R-rest and run main sweeps.
    // bool _rr_graph_warm_start = false;
    // size_t _rr_graph_warm_start_div = 10;
    // size_t _rr_graph_warm_start_sweeps = 1;
    double _hoeffding_delta = 0.05;
    double _hoeffding_margin_scale = 1.0;
    double _hoeffding_margin = 0.1;

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
            else if (!param.compare("-eps")) { _eps = stod(value); _eps_set = true; }
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
            else if (!param.compare("-eval_xi"))
            {
                if (value == "half")
                {
                    _eval_xi_half_set = true;
                    _eval_xi_const_set = true;
                }
                else if (value == "unirand")
                {
                    _eval_xi_unirand_set = true;
                }
                else if (value == "ones" || value == "1")
                {
                    _eval_xi_const = 1.0;
                    _eval_xi_const_set = true;
                }
                else if (value == "lo")
                {
                    _eval_xi_const = _xi_lo;
                    _eval_xi_const_set = true;
                }
                else
                {
                    _eval_xi_const = stod(value);
                    _eval_xi_const_set = true;
                }
            }
            else if (!param.compare("-eval_mc")) _eval_mc = static_cast<uint32_t>(stoul(value));
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
            else if (!param.compare("-q_mode")) _qModeStr = value;
            else if (!param.compare("-q_exp_scale")) _qExpScale = stod(value);
            else if (!param.compare("-q_lo")) _qUnirandLo = stod(value);
            else if (!param.compare("-q_hi")) _qUnirandHi = stod(value);
            else if (!param.compare("-nodehyper_suffix")) _nodehyperSuffix = value;
            else if (!param.compare("-nodehyper_base_suffix")) _nodehyperBaseSuffix = value;
            else if (!param.compare("-hyp_output_suffix")) _hypOutputSuffix = value;
            else if (!param.compare("-hyp_profile")) _hypProfilePath = value;
            else if (!param.compare("-tau_mode")) _tauModeStr = value;
            else if (!param.compare("-tau_lo")) _tauLo = stod(value);
            else if (!param.compare("-tau_hi")) _tauHi = stod(value);
            else if (!param.compare("-tau_jitter")) _tauJitter = stod(value);
            else if (!param.compare("-tau_q_var")) _tauQVar = stod(value);
            else if (!param.compare("-xi_lo")) _xi_lo = stod(value);
            else if (!param.compare("-xi_init"))
            {
                if (value == "lo" || value == "xi_lo")
                    _xi_init_lo = true;
                else if (value == "1" || value == "ones")
                    _xi_init_lo = false;
            }
            else if (!param.compare("-frim_rr"))
            {
                _frim_rr = static_cast<size_t>(stoull(value));
                _frim_rr_fixed = true;
            }
            else if (!param.compare("-frim_mc")) _frim_mc = static_cast<size_t>(stoull(value));
            else if (!param.compare("-frim_sweeps")) _frim_max_sweeps = stoi(value);
            else if (!param.compare("-rr_graph_gate_sweep_index"))
                _rr_graph_gate_sweep_index = (value == "1" || value == "true");
            else if (!param.compare("-rr_graph_keep_no_hit"))
                _rr_graph_store_hit_only = !(value == "1" || value == "true");
            else if (!param.compare("-rr_graph_store_hit_only"))
                _rr_graph_store_hit_only = (value == "1" || value == "true");
            else if (!param.compare("-rr_root_stat_build"))
                _rr_root_stat_build = (value == "1" || value == "true");
#if 0  // warm-start disabled
            else if (!param.compare("-rr_graph_warm_start"))
                _rr_graph_warm_start = (value == "1" || value == "true");
            else if (!param.compare("-rr_graph_warm_start_div"))
                _rr_graph_warm_start_div = static_cast<size_t>(stoull(value));
            else if (!param.compare("-rr_graph_warm_start_sweeps"))
                _rr_graph_warm_start_sweeps = static_cast<size_t>(stoull(value));
#endif
            else if (!param.compare("-hoeffding_delta")) _hoeffding_delta = stod(value);
            else if (!param.compare("-hoeffding_margin_scale"))
                _hoeffding_margin_scale = stod(value);
            else if (!param.compare("-hoeffding_margin")) _hoeffding_margin = stod(value);
        }

        if (_wcVar <= 0)
        {
            //wrong input
            _wcVar = 1.0;
        }

        decode_wrr_sample_mode();
        decode_lam_mode();
        decode_q_mode();
        decode_tau_mode();

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

    std::string pruneMarginSuffix() const
    {
        std::ostringstream oss;
        if (std::abs(_hoeffding_delta - 0.05) >= 1e-12)
            oss << "_d" << _hoeffding_delta;
        if (_hoeffding_margin >= 0.0)
            oss << "_m" << _hoeffding_margin;
        else if (std::abs(_hoeffding_margin_scale - 1.0) >= 1e-12)
            oss << "_ms" << _hoeffding_margin_scale;
        return oss.str();
    }

    std::string frimSweepSuffix() const
    {
        return "_S" + std::to_string(_frim_max_sweeps);
    }

    /// Filename R tag: R0 when -eps was set; otherwise actual sample count.
    std::string frimFilenameRSuffix() const
    {
        if (_eps_set)
            return "_R0";
        return "_R" + std::to_string(_frim_rr);
    }

    /// Filename eps tag when -eps was explicitly passed.
    std::string frimFilenameEpsSuffix() const
    {
        if (!_eps_set)
            return "";
        std::ostringstream oss;
        oss << _eps;
        return "_eps" + oss.str();
    }

    /// RR-graph sample count: fixed -frim_rr, or
    /// ceil((8+2*eps)*T*(log(1/delta)+n*log(2)+log(2))/(eps^2*n))
    /// where T = tau_sum (falls back to n when tau_sum <= 0),
    /// delta = (user delta or 1/n) / 2.
    size_t resolveFrimRRGraphNumSamples(size_t num_v, double tau_sum = -1.0) const
    {
        if (_frim_rr_fixed)
            return _frim_rr;
        if (_eps <= 0.0 || num_v == 0)
            return std::max<size_t>(1, _frim_rr);
        const double n = static_cast<double>(num_v);
        const double T = (tau_sum > 0.0) ? tau_sum : n;
        const double base_delta = (_delta > 0.0) ? _delta : (1.0 / n);
        const double delta = base_delta / 2.0;
        const double eps2 = _eps * _eps;
        const double log_term =
            std::log(1.0 / delta) + n * std::log(2.0) + std::log(2.0);
        const double num_rr =
            (8.0 + 2.0 * _eps) * T * log_term / (eps2 * n);
        return std::max<size_t>(1, static_cast<size_t>(std::ceil(num_rr)));
    }

    void build_frim_outfilename(FuncType func)
    {
        const std::string rTag = frimFilenameRSuffix();
        const std::string epsTag = frimFilenameEpsSuffix();
        if (func == FRIM_RR)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_rr_" + _probDistStr + rTag + epsTag;
        }
        else if (func == FRIM_RR_NAIVE)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_rr_naive_" + _probDistStr + rTag + epsTag
                + frimSweepSuffix();
        }
        else if (func == FRIM_RR_GRAPH)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_rr_graph_" + _probDistStr + rTag + epsTag
                + frimSweepSuffix();
        }
        else if (func == FRIM_RR_GRAPH_PRUNE)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_rr_graph_prune_" + _probDistStr + rTag + epsTag
                + frimSweepSuffix() + pruneMarginSuffix();
        }
        else if (func == FRIM_SUBSIM)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_subsim_" + _probDistStr + epsTag;
        }
        else if (func == FRIM_MC_CRN)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_mc_crn_" + _probDistStr + "_mc" + std::to_string(_frim_mc)
                + frimSweepSuffix();
        }
        else if (func == FRIM_MC_NAIVE)
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed)
                + "_frim_mc_naive_" + _probDistStr + "_mc" + std::to_string(_frim_mc)
                + frimSweepSuffix();
        }
        else
        {
            _outFileName = _graphname + "_" + std::to_string(_rand_seed) + "_frim_" + _probDistStr;
        }
    }

    void build_eval_xi_half_outfilename()
    {
        _outFileName = _graphname + "_" + std::to_string(_rand_seed)
            + "_xi_half_" + _probDistStr;
    }

    void build_eval_xi_unirand_outfilename()
    {
        _outFileName = _graphname + "_" + std::to_string(_rand_seed)
            + "_xi_unirand_" + _probDistStr;
    }

    void build_eval_xi_const_outfilename(double xi_val)
    {
        std::ostringstream oss;
        oss << xi_val;
        _outFileName = _graphname + "_" + std::to_string(_rand_seed)
            + "_xi" + oss.str() + "_" + _probDistStr;
    }

    FuncType decode_frim_eval_method() const
    {
        if (_method == "frim_rr" || _method == "FRIM_RR")
            return FRIM_RR;
        if (_method == "frim_rr_naive" || _method == "FRIM_RR_NAIVE")
            return FRIM_RR_NAIVE;
        if (_method == "frim_rr_graph" || _method == "FRIM_RR_GRAPH"
            || _method == "frim_rr_crn" || _method == "FRIM_RR_CRN")
            return FRIM_RR_GRAPH;
        if (_method == "frim_rr_graph_prune" || _method == "FRIM_RR_GRAPH_PRUNE"
            || _method == "frim_rr_crn_prune" || _method == "FRIM_RR_CRN_PRUNE")
            return FRIM_RR_GRAPH_PRUNE;
        if (_method == "frim_subsim" || _method == "FRIM_SUBSIM")
            return FRIM_SUBSIM;
        if (_method == "frim_mc_crn" || _method == "FRIM_MC_CRN"
            || _method == "frim_mc" || _method == "FRIM_MC")
            return FRIM_MC_CRN;
        if (_method == "frim_mc_naive" || _method == "FRIM_MC_NAIVE")
            return FRIM_MC_NAIVE;
        return FUNC_ERROR;
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
        else if (_funcStr == "hyp_derive" || _funcStr == "HYP_DERIVE")
        {
            _func = HYP_DERIVE;
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
        else if (_funcStr == "frim_rr" || _funcStr == "FRIM_RR")
        {
            _func = FRIM_RR;
        }
        else if (_funcStr == "frim_rr_naive" || _funcStr == "FRIM_RR_NAIVE")
        {
            _func = FRIM_RR_NAIVE;
        }
        else if (_funcStr == "frim_rr_graph" || _funcStr == "FRIM_RR_GRAPH"
                 || _funcStr == "frim_rr_crn" || _funcStr == "FRIM_RR_CRN")
        {
            _func = FRIM_RR_GRAPH;
        }
        else if (_funcStr == "frim_rr_graph_prune" || _funcStr == "FRIM_RR_GRAPH_PRUNE"
                 || _funcStr == "frim_rr_crn_prune" || _funcStr == "FRIM_RR_CRN_PRUNE")
        {
            _func = FRIM_RR_GRAPH_PRUNE;
        }
        else if (_funcStr == "frim_rr_root_stat" || _funcStr == "FRIM_RR_GRAPH_ROOT_STAT")
        {
            _func = FRIM_RR_GRAPH_ROOT_STAT;
        }
        else if (_funcStr == "frim_subsim" || _funcStr == "FRIM_SUBSIM")
        {
            _func = FRIM_SUBSIM;
        }
        else if (_funcStr == "frim_prune" || _funcStr == "FRIM_PRUNE")
        {
            _func = FRIM_PRUNE;
        }
        else if (_funcStr == "frim_prune_lo" || _funcStr == "FRIM_PRUNE_LO")
        {
            _func = FRIM_PRUNE_LO;
        }
        else if (_funcStr == "frim_prune_both" || _funcStr == "FRIM_PRUNE_BOTH")
        {
            _func = FRIM_PRUNE_BOTH;
        }
        else if (_funcStr == "frim_mc_crn" || _funcStr == "FRIM_MC_CRN"
                 || _funcStr == "frim_mc" || _funcStr == "FRIM_MC")
        {
            _func = FRIM_MC_CRN;
        }
        else if (_funcStr == "frim_mc_naive" || _funcStr == "FRIM_MC_NAIVE")
        {
            _func = FRIM_MC_NAIVE;
        }
        else if (_funcStr == "frim_outname" || _funcStr == "FRIM_OUTNAME")
        {
            _func = FRIM_OUTNAME;
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
        if (_lamModeStr == "uniform" || _lamModeStr == "unirand" || _lamModeStr == "two_tier"
            || _lamModeStr == "exponential" || _lamModeStr == "exp")
        {
            _lamModeError = false;
        }
        else
        {
            _lamModeError = true;
        }
    }

    void decode_q_mode()
    {
        if (_qModeStr == "active_inactive" || _qModeStr == "normal"
            || _qModeStr == "exponential" || _qModeStr == "exp"
            || _qModeStr == "unirand" || _qModeStr == "uniform")
        {
            _qModeError = false;
        }
        else
        {
            _qModeError = true;
        }
    }

    void decode_tau_mode()
    {
        if (_tauModeStr == "exponential" || _tauModeStr == "ranked"
            || _tauModeStr == "soft" || _tauModeStr == "exponential_soft"
            || _tauModeStr == "q_normal" || _tauModeStr == "q_plus_normal"
            || _tauModeStr == "uniform"
            || _tauModeStr == "exp_random" || _tauModeStr == "exponential_random")
        {
            _tauModeError = false;
        }
        else
        {
            _tauModeError = true;
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
