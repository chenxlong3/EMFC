#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include "commonStruct.h"

bool greater_first(const std::pair<uint32_t, float> x, const std::pair<uint32_t, float> y)
{
    if (x.second > y.second) return true;
    else if (x.second < y.second) return false;
    else if (x.first < y.first) return true;
    else return false;
}

class GraphBase
{
public:
    enum class LamAssignMode { UNIFORM, UNIRAND, TWO_TIER, EXPONENTIAL };
    enum class QAssignMode { ACTIVE_INACTIVE_NORMAL, EXPONENTIAL, UNIRAND };
    enum class TauAssignMode { UNIFORM_0_1, EXPONENTIAL_NORM, EXPONENTIAL_RANDOM, Q_PLUS_NORMAL };

    /// Configurable assignment modes for per-node hyperparameters (q, tau, lam).
    struct NodeHyperParamsConfig
    {
        LamAssignMode lam_mode = LamAssignMode::UNIRAND;
        double lam_uniform_value = 0.3;
        double lam_high_ratio = 0.3;
        double lam_high_lo = 0.6;
        double lam_high_hi = 0.8;
        double lam_low_lo = 0.1;
        double lam_low_hi = 0.3;

        QAssignMode q_mode = QAssignMode::EXPONENTIAL;
        double q_exponential_scale = 1.0;
        double q_unirand_lo = 0.0;
        double q_unirand_hi = 1.0;
        double active_user_ratio = 0.2;
        double active_q_mean = 0.8;
        double active_q_var = 1.0;
        double inactive_q_mean = 0.2;
        double inactive_q_var = 1.0;

        TauAssignMode tau_mode = TauAssignMode::Q_PLUS_NORMAL;
        double tau_exponential_scale = 1.0;
        double tau_lo = 1.0;
        double tau_hi = 5.0;
        /// Add Uniform(-tau_jitter, +tau_jitter) noise to q before rank-matching tau.
        /// 0 = strict rank; typical soft value ~0.15–0.25 when q in [0,1].
        double tau_jitter = 0.0;

        /// Noise variance for tau_mode=q_normal: tau = clamp((q + N(0, tau_q_var)) * tau_hi, ...).
        double tau_q_var = 1.0;

        double lam_unirand_lo = 0.1;
        double lam_unirand_hi = 0.5;
        double lam_exponential_scale = 1.0;
    };

    static const NodeHyperParamsConfig& defaultNodeHyperParamsConfig()
    {
        static const NodeHyperParamsConfig config{};
        return config;
    }

    static NodeHyperParamsConfig makeNodeHyperParamsConfig(
        const std::string& lam_mode_str,
        double lam_uniform_value,
        double lam_unirand_lo,
        double lam_unirand_hi,
        double active_user_ratio,
        double active_q_mean,
        double active_q_var,
        double inactive_q_mean,
        double inactive_q_var,
        const std::string& q_mode_str = "exponential",
        double q_exp_scale = 1.0,
        double q_unirand_lo = 0.0,
        double q_unirand_hi = 1.0,
        const std::string& tau_mode_str = "q_normal",
        double tau_lo = 1.0,
        double tau_hi = 5.0,
        double tau_jitter = -1.0,
        double tau_q_var = 1.0)
    {
        NodeHyperParamsConfig config;
        config.lam_mode = LamAssignMode::UNIRAND;
        if (lam_mode_str == "uniform")
            config.lam_mode = LamAssignMode::UNIFORM;
        else if (lam_mode_str == "two_tier")
            config.lam_mode = LamAssignMode::TWO_TIER;
        else if (lam_mode_str == "exponential" || lam_mode_str == "exp")
            config.lam_mode = LamAssignMode::EXPONENTIAL;
        config.lam_uniform_value = lam_uniform_value;
        config.lam_unirand_lo = lam_unirand_lo;
        config.lam_unirand_hi = lam_unirand_hi;
        config.active_user_ratio = active_user_ratio;
        config.active_q_mean = active_q_mean;
        config.active_q_var = active_q_var;
        config.inactive_q_mean = inactive_q_mean;
        config.inactive_q_var = inactive_q_var;

        config.q_exponential_scale = q_exp_scale;
        config.q_unirand_lo = q_unirand_lo;
        config.q_unirand_hi = q_unirand_hi;
        if (q_mode_str == "exponential" || q_mode_str == "exp")
            config.q_mode = QAssignMode::EXPONENTIAL;
        else if (q_mode_str == "unirand" || q_mode_str == "uniform")
            config.q_mode = QAssignMode::UNIRAND;
        else
            config.q_mode = QAssignMode::ACTIVE_INACTIVE_NORMAL;

        config.tau_lo = tau_lo;
        config.tau_hi = tau_hi;
        config.tau_q_var = (tau_q_var > 0.0) ? tau_q_var : 1.0;
        if (tau_mode_str == "uniform")
        {
            config.tau_mode = TauAssignMode::UNIFORM_0_1;
            config.tau_jitter = 0.0;
        }
        else if (tau_mode_str == "exp_random" || tau_mode_str == "exponential_random")
        {
            config.tau_mode = TauAssignMode::EXPONENTIAL_RANDOM;
            config.tau_jitter = 0.0;
        }
        else if (tau_mode_str == "soft" || tau_mode_str == "exponential_soft")
        {
            config.tau_mode = TauAssignMode::EXPONENTIAL_NORM;
            config.tau_jitter = (tau_jitter >= 0.0) ? tau_jitter : 0.2;
        }
        else if (tau_mode_str == "q_normal" || tau_mode_str == "q_plus_normal")
        {
            config.tau_mode = TauAssignMode::Q_PLUS_NORMAL;
            config.tau_jitter = 0.0;
        }
        else
        {
            config.tau_mode = TauAssignMode::EXPONENTIAL_NORM;
            config.tau_jitter = (tau_jitter >= 0.0) ? tau_jitter : 0.0;
        }
        return config;
    }

    struct NodeHyperParams
    {
        std::vector<double> q;
        std::vector<double> tau;
        std::vector<double> lam;
    };

    /// Load baseline nodehyper.vec, override selected fields, save to a suffix file.
    struct HypDeriveSpec
    {
        std::string base_suffix;
        std::string output_suffix;
        bool keep_tau = true;
        bool keep_lam = true;
        bool override_q = false;
        bool override_lam = false;
        bool override_tau = false;
        int rand_seed = -1;
        NodeHyperParamsConfig q_config;
        NodeHyperParamsConfig lam_config;
        NodeHyperParamsConfig tau_config;
    };

    static std::string trimHypProfileValue(const std::string& s)
    {
        const size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
            return "";
        const size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

    static bool parseHypProfileBool(const std::string& value, bool default_value)
    {
        const std::string v = trimHypProfileValue(value);
        if (v.empty())
            return default_value;
        if (v == "1" || v == "true" || v == "yes" || v == "on")
            return true;
        if (v == "0" || v == "false" || v == "no" || v == "off")
            return false;
        return default_value;
    }

    static bool loadHypDeriveProfile(const std::string& profile_path, HypDeriveSpec& spec)
    {
        std::ifstream infile(profile_path);
        if (!infile.is_open())
        {
            std::cout << "Cannot open hyp profile: " << profile_path << std::endl;
            return false;
        }

        std::string q_mode = "keep";
        std::string lam_mode = "keep";
        std::string tau_mode = "keep";
        double q_lo = 0.0;
        double q_hi = 1.0;
        double q_exp_scale = 1.0;
        double lam_lo = 0.1;
        double lam_hi = 0.5;
        double lam_value = 0.3;
        double lam_exp_scale = 1.0;
        double tau_lo = 1.0;
        double tau_hi = 5.0;
        double tau_jitter = -1.0;
        double tau_q_var = 1.0;

        std::string line;
        while (std::getline(infile, line))
        {
            line = trimHypProfileValue(line);
            if (line.empty() || line[0] == '#')
                continue;
            const size_t eq = line.find('=');
            if (eq == std::string::npos)
                continue;
            const std::string key = trimHypProfileValue(line.substr(0, eq));
            const std::string value = trimHypProfileValue(line.substr(eq + 1));
            if (key == "base_suffix")
                spec.base_suffix = value;
            else if (key == "output_suffix")
                spec.output_suffix = value;
            else if (key == "keep_tau")
                spec.keep_tau = parseHypProfileBool(value, spec.keep_tau);
            else if (key == "keep_lam")
                spec.keep_lam = parseHypProfileBool(value, spec.keep_lam);
            else if (key == "rand_seed")
                spec.rand_seed = std::stoi(value);
            else if (key == "q_mode")
                q_mode = value;
            else if (key == "q_lo")
                q_lo = std::stod(value);
            else if (key == "q_hi")
                q_hi = std::stod(value);
            else if (key == "q_exp_scale")
                q_exp_scale = std::stod(value);
            else if (key == "lam_mode")
                lam_mode = value;
            else if (key == "lam_lo")
                lam_lo = std::stod(value);
            else if (key == "lam_hi")
                lam_hi = std::stod(value);
            else if (key == "lam_value")
                lam_value = std::stod(value);
            else if (key == "lam_exp_scale")
                lam_exp_scale = std::stod(value);
            else if (key == "tau_mode")
                tau_mode = value;
            else if (key == "tau_lo")
                tau_lo = std::stod(value);
            else if (key == "tau_hi")
                tau_hi = std::stod(value);
            else if (key == "tau_jitter")
                tau_jitter = std::stod(value);
            else if (key == "tau_q_var")
                tau_q_var = std::stod(value);
        }

        if (spec.output_suffix.empty())
        {
            std::cout << "hyp profile missing output_suffix: " << profile_path << std::endl;
            return false;
        }

        if (q_mode != "keep")
        {
            spec.override_q = true;
            spec.q_config = makeNodeHyperParamsConfig(
                "unirand", lam_value, lam_lo, lam_hi,
                0.2, 0.8, 1.0, 0.2, 1.0,
                q_mode, q_exp_scale, q_lo, q_hi,
                "q_normal", tau_lo, tau_hi, tau_jitter, tau_q_var);
        }
        if (!spec.keep_lam && lam_mode != "keep")
        {
            spec.override_lam = true;
            spec.lam_config = makeNodeHyperParamsConfig(
                lam_mode, lam_value, lam_lo, lam_hi,
                0.2, 0.8, 1.0, 0.2, 1.0,
                "exponential", q_exp_scale, q_lo, q_hi,
                "q_normal", tau_lo, tau_hi, tau_jitter, tau_q_var);
            spec.lam_config.lam_exponential_scale = lam_exp_scale;
        }
        if (!spec.keep_tau && tau_mode != "keep")
        {
            spec.override_tau = true;
            spec.tau_config = makeNodeHyperParamsConfig(
                "unirand", lam_value, lam_lo, lam_hi,
                0.2, 0.8, 1.0, 0.2, 1.0,
                "exponential", q_exp_scale, q_lo, q_hi,
                tau_mode, tau_lo, tau_hi, tau_jitter, tau_q_var);
        }
        return true;
    }

    static bool deriveAndSaveNodeHyperParams(
        const std::string& graphName,
        const HypDeriveSpec& spec)
    {
        NodeHyperParams params;
        if (!TIO::readGraphNodeHyperParamsFileWithSuffix(
                graphName, spec.base_suffix, params.q, params.tau, params.lam))
        {
            std::cout << "Cannot load base node hyper params for " << graphName;
            if (!spec.base_suffix.empty())
                std::cout << " (suffix=" << spec.base_suffix << ")";
            std::cout << std::endl;
            return false;
        }

        NodeHyperParams scratch = params;
        if (spec.override_q)
        {
            assignQ(spec.q_config, scratch);
            params.q = scratch.q;
        }
        if (spec.override_lam)
        {
            assignLam(spec.lam_config, scratch);
            params.lam = scratch.lam;
        }
        if (spec.override_tau)
        {
            assignTau(spec.tau_config, scratch);
            params.tau = scratch.tau;
        }

        TIO::SaveGraphNodeHyperParamsWithSuffix(
            graphName, spec.output_suffix, params.q, params.tau, params.lam);

        double q_sum = 0.0;
        double tau_sum = 0.0;
        double lam_sum = 0.0;
        for (size_t i = 0; i < params.q.size(); i++)
        {
            q_sum += params.q[i];
            tau_sum += params.tau[i];
            lam_sum += params.lam[i];
        }
        std::cout << "Derived node hyper params -> "
                  << graphName << ".nodehyper." << spec.output_suffix << ".vec"
                  << " (n=" << params.q.size()
                  << ", sum q=" << q_sum
                  << ", sum tau=" << tau_sum
                  << ", sum lam=" << lam_sum << ")" << std::endl;
        return true;
    }

    /// Directed graph: graph[u] = out-edges from u. avg_deg = |E|/|V|; four Pearson r on each directed edge u->v.
    struct DegreeAssortativityStats
    {
        uint64_t num_dir_edges;              ///< number of directed edges |E|
        double avg_deg;                      ///< |E| / |V| (same as mean in- or out-degree)
        double assortativity_out_out;        ///< corr(k_out(u), k_out(v))
        double assortativity_out_in;         ///< corr(k_out(u), k_in(v))
        double assortativity_in_out;         ///< corr(k_in(u), k_out(v))
        double assortativity_in_in;          ///< corr(k_in(u), k_in(v))
    };

    /// Compute avg_deg and four directed degree assortativities (Pearson r over directed edges).
    static DegreeAssortativityStats computeDegreeAssortativityStats(const Graph& graph)
    {
        DegreeAssortativityStats s{};
        const size_t n = graph.size();
        if (n == 0)
            return s;

        std::vector<size_t> in_deg(n, 0);
        size_t m = 0;
        for (size_t u = 0; u < n; u++)
        {
            m += graph[u].size();
            for (const auto& e : graph[u])
                in_deg[e.first]++;
        }

        s.num_dir_edges = static_cast<uint64_t>(m);
        s.avg_deg = static_cast<double>(m) / static_cast<double>(n);
        if (m == 0)
            return s;

        std::vector<size_t> out_deg(n);
        for (size_t u = 0; u < n; u++)
            out_deg[u] = graph[u].size();

        double sx_oo = 0.0, sy_oo = 0.0, sxx_oo = 0.0, syy_oo = 0.0, sxy_oo = 0.0;
        double sx_oi = 0.0, sy_oi = 0.0, sxx_oi = 0.0, syy_oi = 0.0, sxy_oi = 0.0;
        double sx_io = 0.0, sy_io = 0.0, sxx_io = 0.0, syy_io = 0.0, sxy_io = 0.0;
        double sx_ii = 0.0, sy_ii = 0.0, sxx_ii = 0.0, syy_ii = 0.0, sxy_ii = 0.0;

        for (size_t u = 0; u < n; u++)
        {
            const double ou = static_cast<double>(out_deg[u]);
            const double iu = static_cast<double>(in_deg[u]);
            for (const auto& e : graph[u])
            {
                const size_t v = e.first;
                const double ov = static_cast<double>(out_deg[v]);
                const double iv = static_cast<double>(in_deg[v]);

                sx_oo += ou;
                sy_oo += ov;
                sxx_oo += ou * ou;
                syy_oo += ov * ov;
                sxy_oo += ou * ov;

                sx_oi += ou;
                sy_oi += iv;
                sxx_oi += ou * ou;
                syy_oi += iv * iv;
                sxy_oi += ou * iv;

                sx_io += iu;
                sy_io += ov;
                sxx_io += iu * iu;
                syy_io += ov * ov;
                sxy_io += iu * ov;

                sx_ii += iu;
                sy_ii += iv;
                sxx_ii += iu * iu;
                syy_ii += iv * iv;
                sxy_ii += iu * iv;
            }
        }

        const double md = static_cast<double>(m);
        s.assortativity_out_out = pearsonFromSums(sx_oo, sy_oo, sxx_oo, syy_oo, sxy_oo, md);
        s.assortativity_out_in = pearsonFromSums(sx_oi, sy_oi, sxx_oi, syy_oi, sxy_oi, md);
        s.assortativity_in_out = pearsonFromSums(sx_io, sy_io, sxx_io, syy_io, sxy_io, md);
        s.assortativity_in_in = pearsonFromSums(sx_ii, sy_ii, sxx_ii, syy_ii, sxy_ii, md);

        return s;
    }

    /// Append avg_deg and four assortativity lines (key: value) to graphInfo/{dataset_name}.attr
    static void appendDegreeStatsToAttrFile(
        const std::string& dataset_name,
        const DegreeAssortativityStats& stats,
        const std::string& graph_info_dir = "graphInfo")
    {
        const std::string path = graph_info_dir + "/" + dataset_name + ".attr";
        std::ofstream ofs(path, std::ios::app);
        if (!ofs.is_open())
        {
            std::cerr << "appendDegreeStatsToAttrFile: cannot open " << path << std::endl;
            return;
        }
        ofs << "num dir edges: " << stats.num_dir_edges << "\n";
        ofs << std::fixed << std::setprecision(12);
        ofs << "avg_deg: " << stats.avg_deg << "\n";
        ofs << "assortativity_out_out: " << stats.assortativity_out_out << "\n";
        ofs << "assortativity_out_in: " << stats.assortativity_out_in << "\n";
        ofs << "assortativity_in_out: " << stats.assortativity_in_out << "\n";
        ofs << "assortativity_in_in: " << stats.assortativity_in_in << "\n";
    }

    /// Compute stats and append to graphInfo/{dataset_name}.attr (same keys as appendDegreeStatsToAttrFile).
    static void computeAndAppendDegreeStatsToAttr(
        const Graph& graph,
        const std::string& dataset_name,
        const std::string& graph_info_dir = "graphInfo")
    {
        const DegreeAssortativityStats st = computeDegreeAssortativityStats(graph);
        appendDegreeStatsToAttrFile(dataset_name, st, graph_info_dir);
    }

    static double pearsonFromSums(double sx, double sy, double sxx, double syy, double sxy, double md)
    {
        if (md <= 0.0)
            return 0.0;
        const double mx = sx / md;
        const double my = sy / md;
        const double cov = sxy / md - mx * my;
        const double vx = sxx / md - mx * mx;
        const double vy = syy / md - my * my;
        if (vx <= 0.0 || vy <= 0.0)
            return 0.0;
        return cov / (std::sqrt(vx) * std::sqrt(vy));
    }

    static double clampUnitInterval(double x)
    {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    static double sampleNormal(double mean, double variance)
    {
        static bool has_spare = false;
        static double spare = 0.0;
        const double stddev = std::sqrt(variance);

        if (has_spare)
        {
            has_spare = false;
            return mean + stddev * spare;
        }

        has_spare = true;
        const double u1 = dsfmt_gv_genrand_open_open();
        const double u2 = dsfmt_gv_genrand_open_open();
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * M_PI * u2;
        spare = radius * std::sin(theta);
        return mean + stddev * radius * std::cos(theta);
    }

    static void assignLam(const NodeHyperParamsConfig& config, NodeHyperParams& params)
    {
        const size_t n = params.lam.size();
        switch (config.lam_mode)
        {
        case LamAssignMode::UNIFORM:
            std::fill(params.lam.begin(), params.lam.end(), config.lam_uniform_value);
            break;
        case LamAssignMode::UNIRAND:
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_close_open();
                params.lam[i] = config.lam_unirand_lo + u * (config.lam_unirand_hi - config.lam_unirand_lo);
            }
            break;
        case LamAssignMode::TWO_TIER:
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_close_open();
                if (dsfmt_gv_genrand_open_close() < config.lam_high_ratio)
                    params.lam[i] = config.lam_high_lo + u * (config.lam_high_hi - config.lam_high_lo);
                else
                    params.lam[i] = config.lam_low_lo + u * (config.lam_low_hi - config.lam_low_lo);
            }
            break;
        case LamAssignMode::EXPONENTIAL:
        {
            // Truncated exponential on [0,1]: f(lam) propto exp(-lam / scale).
            const double scale = (config.lam_exponential_scale > 0.0)
                ? config.lam_exponential_scale
                : 1.0;
            const double trunc_mass = 1.0 - std::exp(-1.0 / scale);
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_open_open();
                params.lam[i] = -scale * std::log(1.0 - u * trunc_mass);
            }
            break;
        }
        }
    }

    static double clampRange(double x, double lo, double hi)
    {
        if (x < lo) return lo;
        if (x > hi) return hi;
        return x;
    }

    static void assignTau(const NodeHyperParamsConfig& config, NodeHyperParams& params)
    {
        const size_t n = params.tau.size();
        const double tau_span = config.tau_hi - config.tau_lo;
        switch (config.tau_mode)
        {
        case TauAssignMode::Q_PLUS_NORMAL:
            for (size_t i = 0; i < n; i++)
            {
                const double z = sampleNormal(0.0, config.tau_q_var);
                params.tau[i] = clampRange(
                    (params.q[i] + z) * config.tau_hi,
                    config.tau_lo,
                    config.tau_hi);
            }
            break;
        case TauAssignMode::UNIFORM_0_1:
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_close_open();
                params.tau[i] = config.tau_lo + tau_span * u;
            }
            break;
        case TauAssignMode::EXPONENTIAL_NORM:
        case TauAssignMode::EXPONENTIAL_RANDOM:
        {
            std::vector<double> tau_samples(n);
            double max_tau = 0.0;
            for (size_t i = 0; i < n; i++)
            {
                tau_samples[i] = -config.tau_exponential_scale * std::log(1.0 - dsfmt_gv_genrand_open_open());
                max_tau = std::max(max_tau, tau_samples[i]);
            }
            if (max_tau > 0.0)
            {
                for (size_t i = 0; i < n; i++)
                    tau_samples[i] /= max_tau;
            }
            else
            {
                std::fill(tau_samples.begin(), tau_samples.end(), 0.0);
            }

            for (size_t i = 0; i < n; i++)
                tau_samples[i] = config.tau_lo + tau_span * tau_samples[i];

            std::sort(tau_samples.begin(), tau_samples.end(), std::greater<double>());

            std::vector<uint32_t> order(n);
            std::iota(order.begin(), order.end(), 0U);
            if (config.tau_mode == TauAssignMode::EXPONENTIAL_RANDOM)
            {
                for (size_t i = n - 1; i > 0; i--)
                {
                    const size_t j =
                        dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
                    std::swap(order[i], order[j]);
                }
            }
            else if (config.tau_jitter > 0.0)
            {
                std::vector<double> keys(n);
                for (size_t i = 0; i < n; i++)
                {
                    keys[i] = params.q[i]
                        + config.tau_jitter
                              * (2.0 * dsfmt_gv_genrand_close_open() - 1.0);
                }
                std::sort(
                    order.begin(),
                    order.end(),
                    [&](uint32_t a, uint32_t b) { return keys[a] > keys[b]; });
            }
            else
            {
                std::sort(
                    order.begin(),
                    order.end(),
                    [&](uint32_t a, uint32_t b) { return params.q[a] > params.q[b]; });
            }

            for (size_t rank = 0; rank < n; rank++)
                params.tau[order[rank]] = tau_samples[rank];
            break;
        }
        }
    }

    static void assignQ(const NodeHyperParamsConfig& config, NodeHyperParams& params)
    {
        const size_t n = params.q.size();
        if (n == 0)
            return;

        switch (config.q_mode)
        {
        case QAssignMode::ACTIVE_INACTIVE_NORMAL:
        {
            std::vector<uint32_t> node_ids(n);
            std::iota(node_ids.begin(), node_ids.end(), 0U);
            for (size_t i = n - 1; i > 0; i--)
            {
                const size_t j = dsfmt_gv_genrand_uint32_range(static_cast<uint32_t>(i + 1));
                std::swap(node_ids[i], node_ids[j]);
            }

            const size_t num_active = std::max<size_t>(
                1,
                static_cast<size_t>(std::ceil(config.active_user_ratio * static_cast<double>(n))));
            std::vector<bool> is_active(n, false);
            for (size_t i = 0; i < num_active; i++)
                is_active[node_ids[i]] = true;

            for (size_t i = 0; i < n; i++)
            {
                if (is_active[i])
                    params.q[i] = clampUnitInterval(sampleNormal(config.active_q_mean, config.active_q_var));
                else
                    params.q[i] = clampUnitInterval(sampleNormal(config.inactive_q_mean, config.inactive_q_var));
            }
            break;
        }
        case QAssignMode::EXPONENTIAL:
        {
            // Truncated exponential on [0,1]: f(q) propto exp(-q / scale).
            const double scale = (config.q_exponential_scale > 0.0)
                ? config.q_exponential_scale
                : 1.0;
            const double trunc_mass = 1.0 - std::exp(-1.0 / scale);
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_open_open();
                params.q[i] = -scale * std::log(1.0 - u * trunc_mass);
            }
            break;
        }
        case QAssignMode::UNIRAND:
            for (size_t i = 0; i < n; i++)
            {
                const double u = dsfmt_gv_genrand_close_open();
                params.q[i] = clampUnitInterval(
                    config.q_unirand_lo + u * (config.q_unirand_hi - config.q_unirand_lo));
            }
            break;
        }
    }

    /// Prepare per-node hyperparameters: q (push prob), tau (node value), lam (discount).
    static NodeHyperParams prepareNodeHyperParams(
        size_t num_nodes,
        const NodeHyperParamsConfig& config = defaultNodeHyperParamsConfig())
    {
        NodeHyperParams params;
        params.q.assign(num_nodes, 0.0);
        params.tau.assign(num_nodes, 0.0);
        params.lam.assign(num_nodes, 0.0);

        assignQ(config, params);
        assignTau(config, params);
        assignLam(config, params);
        return params;
    }

    /// Format the input for future computing, which is much faster for loading. Vector serialization is used.
    static void FormatGraph(
        const std::string filename,
        ProbDist probDist,
        const float sum,
        const float prob,
        const std::string skewType,
        const NodeHyperParamsConfig& hyperConfig = defaultNodeHyperParamsConfig())
    {
        size_t numV, numE;
        uint32_t srcId, dstId;
        float weight = 0.0;
        std::ifstream infile(filename);

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return;
        }

        infile >> numV >> numE;
        Graph vecGRev(numV);
        std::vector<size_t> vecInDeg(numV);

        for (auto i = numE; i--;)
        {
            if (probDist == WEIGHTS)
            {
                infile >> srcId >> dstId >> weight;
            }
            else
            {
                infile >> srcId >> dstId;
            }

            vecGRev[dstId].push_back(Edge(srcId, weight));
        }

        infile.close();

        for (auto idx = 0; idx < numV; idx++)
        {
            vecInDeg[idx] = vecGRev[idx].size();
        }

        if (probDist == WC)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                if (vecGRev[i].size() == 0) continue;

                weight = sum / vecInDeg[i];

                for (size_t j = 0; j < vecGRev[i].size(); j++)
                {
                    vecGRev[i][j].second = weight;
                }
            }
        }
        else if (probDist == UNIFORM)
        {
            // Uniform probability
            for (auto& nbrs : vecGRev)
            {
                for (auto& nbr : nbrs)
                {
                    nbr.second = prob;
                }
            }
        }
        // exponential distribution with lambada equal to its in-degree
        else if (probDist == SKEWED)
        {

            if (skewType == "weibull")
            {
                std::default_random_engine generator(time(NULL));
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;

                for (size_t i = 0; i< vecGRev.size(); i++)
                {

                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // random number from (0, 10)
                        double a = dsfmt_gv_genrand_open_open() * 10;
                        double b = dsfmt_gv_genrand_open_open() * 10;
                        std::weibull_distribution<double> distribution(a,b);
                        auto weight =  distribution(generator);
                        vecGRev[i][j].second = weight;

                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        auto weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value;
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
            else
            {
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;
                for (size_t i = 0; i<vecGRev.size(); i++)
                {
                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // lambda = 1
                        auto weight = -log ( 1.0 - dsfmt_gv_genrand_open_open());
                        vecGRev[i][j].second = weight;
                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        double weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value;
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
        }
        else if (probDist == WEIGHTS)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
            }
        }

        // Build forward graph from finalized reverse graph.
        Graph vecG(numV);
        for (size_t dst = 0; dst < vecGRev.size(); dst++)
        {
            for (const auto& e : vecGRev[dst])
            {
                vecG[e.first].push_back(Edge(static_cast<uint32_t>(dst), e.second));
            }
        }

        std::cout << "probability distribution: " << probDist << std::endl;
        TIO::SaveGraphStruct(filename, vecG, false);  // Save forward graph
        TIO::SaveGraphStruct(filename, vecGRev, true);  // Save reverse graph
        TIO::SaveGraphProbDist(filename, (int)probDist);

        const NodeHyperParams nodeHyperParams = prepareNodeHyperParams(numV, hyperConfig);
        TIO::SaveGraphNodeHyperParams(filename, nodeHyperParams.q, nodeHyperParams.tau, nodeHyperParams.lam);

        std::cout << "The graph is formatted!" << std::endl;
    }

    static NodeHyperParams LoadGraphNodeHyperParams(const std::string& graphName)
    {
        NodeHyperParams params;
        TIO::LoadGraphNodeHyperParams(graphName, params.q, params.tau, params.lam);
        return params;
    }

    /// Load reverse graph via vector deserialization.
    static void LoadGraph(Graph &graph, const std::string graphName)
    {
        TIO::LoadGraphStruct(graphName, graph, true);
        return ;
    }

    /// Load forward graph via vector deserialization.
    static void LoadForwardGraph(Graph &graph, const std::string graphName)
    {
        TIO::LoadGraphStruct(graphName, graph, false);
        return ;
    }

    static int LoadGraphProbDist(const std::string graphName)
    {
        return TIO::LoadGraphProbDist(graphName);
    }

    /// Influence evaluation using Monte Carlo simulation
    static double inf_eval(const Graph& graph, const std::vector<uint32_t>& vecSeed, const CascadeModel model,
                           uint32_t evalsize = 10000,
                           bool printinfo=true)
    {
        Timer EvalTimer;
        if (evalsize == 0) evalsize=10000;
        uint32_t nodeId, currProgress = 0;
        const auto numV = graph.size();
        std::queue<uint32_t> Que;
        std::vector<uint32_t> vecActivated;

        double inf = (double)vecSeed.size();
        bool* activated = (bool *)calloc(numV, sizeof(bool));
        uint32_t* visited = (uint32_t *)calloc(numV, sizeof(uint32_t));
        std::vector<double> vecThr(numV);
        std::vector<double> vecActivateWeight(numV, 0.0);
        for (auto seedId : vecSeed) {
            activated[seedId] = true;
        }
        for (uint32_t i = 0; i < evalsize; i++)
        {
            double cur_round_inf = 0.0;
            if (i * 100 >= evalsize * currProgress)
            {
                const auto evalTime = EvalTimer.get_operation_time();
                if (evalTime > 100)
                    std::cout << "\tMC-Progress at: " << currProgress << "%, " << "time used: " << evalTime << std::endl;
                currProgress += 20;
            }
            for (auto seed : vecSeed)
            {
                Que.push(seed);
            }

            // BFS traversal
            if (model == IC)
            {
                while (!Que.empty())
                {
                    nodeId = Que.front();
                    Que.pop();
                    for (auto& nbr : graph[nodeId])
                    {
                        if (activated[nbr.first]) continue;
                        if (dsfmt_gv_genrand_open_close() <= nbr.second)
                        {
                            activated[nbr.first] = true;
                            vecActivated.push_back(nbr.first);
                            Que.push(nbr.first);
                            cur_round_inf += 1.0;
                        }
                    }
                }
            }
            else if (model == LT)
            {
                while (!Que.empty())
                {
                    nodeId = Que.front();
                    Que.pop();
                    for (auto& nbr : graph[nodeId])
                    {
                        if (activated[nbr.first]) continue;
                        if (visited[nbr.first] < i + 1)
                        {
                            // First time visit this node
                            visited[nbr.first] = i + 1;
                            vecThr[nbr.first] = dsfmt_gv_genrand_open_close();
                            vecActivateWeight[nbr.first] = 0.0;
                        }
                        vecActivateWeight[nbr.first] += nbr.second;
                        if (vecActivateWeight[nbr.first] >= vecThr[nbr.first])
                        {
                            // Activation weight is greater than threshold
                            activated[nbr.first] = true;
                            vecActivated.push_back(nbr.first);
                            Que.push(nbr.first);
                            cur_round_inf += 1.0;
                        }
                    }
                }
            }
            inf += cur_round_inf / evalsize;
            for (auto activatedNode : vecActivated) activated[activatedNode] = false;
            vecActivated.clear();
        }
        free(activated);
        free(visited);
        if (printinfo)
        {
            std::cout << "  >>>MC-Influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << std::endl;
        }

        return inf;
    }

    /// Sample a seed set S₀: include u iff U(0,1) <= q[u] (one draw per node).
    static Nodelist sampleSeedSetFromQ(const std::vector<double>& q)
    {
        Nodelist seeds;
        for (size_t u = 0; u < q.size(); u++)
        {
            if (dsfmt_gv_genrand_open_close() <= q[u])
                seeds.push_back(static_cast<uint32_t>(u));
        }
        return seeds;
    }

    /// MC benefit evaluation: sample seeds from q, propagate with IC/LT,
    /// gate forwarding by xi; benefit is tau if not forwarding, tau*(1-lam) if forwarding.
    static double benefit_inf_eval(
        const Graph& graph,
        const std::vector<double>& q,
        const std::vector<double>& tau,
        const std::vector<double>& lam,
        const CascadeModel model,
        const std::vector<double>& xi = {},
        uint32_t evalsize = 10000,
        bool printinfo = true)
    {
        Timer EvalTimer;
        if (evalsize == 0)
            evalsize = 10000;

        const auto numV = graph.size();
        if (numV == 0 || q.size() != numV || tau.size() != numV || lam.size() != numV)
            return 0.0;

        const bool useXi = (xi.size() == numV);
        uint32_t nodeId, currProgress = 0;
        std::queue<uint32_t> Que;
        std::vector<uint32_t> vecActivated;
        bool* activated = static_cast<bool*>(calloc(numV, sizeof(bool)));
        uint32_t* visited = static_cast<uint32_t*>(calloc(numV, sizeof(uint32_t)));
        std::vector<double> vecThr(numV);
        std::vector<double> vecActivateWeight(numV, 0.0);

        double benefit = 0.0;
        for (uint32_t i = 0; i < evalsize; i++)
        {
            double cur_round_benefit = 0.0;
            if (i * 100 >= evalsize * currProgress)
            {
                const auto evalTime = EvalTimer.get_operation_time();
                if (evalTime > 100)
                    std::cout << "\tBenefit-MC-Progress at: " << currProgress << "%, "
                              << "time used: " << evalTime << std::endl;
                currProgress += 20;
            }

            for (size_t u = 0; u < numV; u++)
            {
                if (dsfmt_gv_genrand_open_close() > q[u])
                    continue;
                activated[u] = true;
                vecActivated.push_back(static_cast<uint32_t>(u));
                Que.push(static_cast<uint32_t>(u));
            }

            if (model == IC)
            {
                while (!Que.empty())
                {
                    nodeId = Que.front();
                    Que.pop();

                    const double xi_u = useXi ? xi[nodeId] : 1.0;
                    if (dsfmt_gv_genrand_open_close() <= xi_u)
                    {
                        cur_round_benefit += tau[nodeId] * (1.0 - lam[nodeId]);
                        for (auto& nbr : graph[nodeId])
                        {
                            if (activated[nbr.first])
                                continue;
                            if (dsfmt_gv_genrand_open_close() <= nbr.second)
                            {
                                activated[nbr.first] = true;
                                vecActivated.push_back(nbr.first);
                                Que.push(nbr.first);
                            }
                        }
                    }
                    else
                    {
                        cur_round_benefit += tau[nodeId];
                    }
                }
            }
            else if (model == LT)
            {
                while (!Que.empty())
                {
                    nodeId = Que.front();
                    Que.pop();

                    const double xi_u = useXi ? xi[nodeId] : 1.0;
                    if (dsfmt_gv_genrand_open_close() <= xi_u)
                    {
                        cur_round_benefit += tau[nodeId] * (1.0 - lam[nodeId]);
                        for (auto& nbr : graph[nodeId])
                        {
                            if (activated[nbr.first])
                                continue;
                            if (visited[nbr.first] < i + 1)
                            {
                                visited[nbr.first] = i + 1;
                                vecThr[nbr.first] = dsfmt_gv_genrand_open_close();
                                vecActivateWeight[nbr.first] = 0.0;
                            }
                            vecActivateWeight[nbr.first] += nbr.second;
                            if (vecActivateWeight[nbr.first] >= vecThr[nbr.first])
                            {
                                activated[nbr.first] = true;
                                vecActivated.push_back(nbr.first);
                                Que.push(nbr.first);
                            }
                        }
                    }
                    else
                    {
                        cur_round_benefit += tau[nodeId];
                    }
                }
            }

            benefit += cur_round_benefit / evalsize;
            for (auto activatedNode : vecActivated)
                activated[activatedNode] = false;
            vecActivated.clear();
        }

        free(activated);
        free(visited);

        if (printinfo)
        {
            std::cout << "  >>>MC-Benefit: " << benefit
                      << ", time used (sec): " << EvalTimer.get_total_time() << std::endl;
        }

        return benefit;
    }
};
