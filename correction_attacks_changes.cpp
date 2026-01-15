// Our modifications were implemented on 2025-10-03.
//
// This file is an improved variant of the original implementation.
// Original-author reference:
// - Repository: https://github.com/Chair-for-Security-Engineering/dilithium-faults
// - File: skippingFault/correction_attacks.cpp
//
// What we changed (relative to the original correction_attacks.cpp by Krahmer et al.):
// 1) Added an ILP backend using Gurobi: solve_ilp_gurobi(...).
// 2) Added two compile-time knobs:
//    - ILP_TIME_LIMIT_MS: ILP solver time limit in milliseconds.
//    - MIN_ROWS_TO_TRY_ILP: minimum number of collected equations before attempting ILP.
// 3) Modified the three recovery functions
//    - recover_s1_skip_comp(...)
//    - recover_s1_skip(...)
//    - recover_s1_skip_shuff(...)
//    by replacing the original NTL::solve()-based linear-system solve with an ILP/Gurobi path.
//    In particular, we attempt ILP as soon as collected >= MIN_ROWS_TO_TRY_ILP, so recovery may
//    succeed before collecting the full N equations.
//
// IMPORTANT include-order note (Gurobi):
// - Keep #include "gurobi_c++.h" BEFORE any extern "C" blocks.
//   Putting a C++ header under extern "C" may cause C linkage to be applied incorrectly,
//   leading to signature mismatches, re-declarations, or link errors.

#include "gurobi_c++.h"
#include <cmath>
#include <limits>
#include <stdexcept>

// ILP solver time limit (milliseconds)
#ifndef ILP_TIME_LIMIT_MS
#define ILP_TIME_LIMIT_MS 60000
#endif

// Heuristic threshold (number of equations) before attempting ILP.
#ifndef MIN_ROWS_TO_TRY_ILP
#define MIN_ROWS_TO_TRY_ILP 193
#endif

bool solve_ilp_gurobi(
    NTL::Vec<NTL::ZZ_p>& solution,
    const NTL::Mat<NTL::ZZ_p>& eq_sys,
    const NTL::Vec<NTL::ZZ_p>& eq_sys_b,
    int m_eq,
    int n_var,
    int eta,
    const std::vector<bool>& exact_rows,
    int time_limit_ms = ILP_TIME_LIMIT_MS,
    int verbose = 0)
{
    using namespace NTL;
    try {
        // Initialize Gurobi.
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_OutputFlag, verbose);
        env.set(GRB_DoubleParam_TimeLimit, (double)time_limit_ms / 1000.0);
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_StringAttr_ModelName, "ILP_Recover_S1");

        // Variables: x_j ∈ [-eta, eta] (integer).
        std::vector<GRBVar> x_vars;
        x_vars.reserve(n_var);
        for (int j = 0; j < n_var; j++) {
            GRBVar x = model.addVar(-eta, eta, 0.0, GRB_INTEGER, "x_" + std::to_string(j));
            x_vars.push_back(x);
        }

        // Add linear constraints: A * x = b over centered representatives in Z.
        for (int i = 0; i < m_eq; i++) {
            if (!exact_rows[i]) continue;
            GRBLinExpr expr = 0.0;
            for (int j = 0; j < n_var; j++) {
                long coeff = conv<long>(rep(eq_sys[i][j]));
                if (coeff != 0) {
                    if (coeff > Q / 2) coeff -= Q;  // center mod
                    expr += coeff * x_vars[j];
                }
            }
            long rhs = conv<long>(rep(eq_sys_b[i]));
            if (rhs > Q / 2) rhs -= Q;
            model.addConstr(expr == rhs, "eq_" + std::to_string(i));
        }

        // Select an objective function.
        GRBLinExpr obj = 0.0;
        for (auto &x : x_vars) obj += x;
        model.setObjective(obj, GRB_MINIMIZE);

        // Solve for solutions.
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            solution.SetLength(n_var);
            for (int j = 0; j < n_var; j++) {
                long val = std::llround(x_vars[j].get(GRB_DoubleAttr_X));
                long modv = ((val % Q) + Q) % Q;
                solution[j] = ZZ_p(modv);
            }
            if (verbose)
                std::cerr << "[ILP-Gurobi] Solved: Obj=" << model.get(GRB_DoubleAttr_ObjVal)
                          << " Time=" << model.get(GRB_DoubleAttr_Runtime) << "s\n";
            return true;
        } else {
            if (verbose)
                std::cerr << "[ILP-Gurobi] Failed: status=" << status << std::endl;
            return false;
        }

    } catch (GRBException &e) {
        if (verbose)
            std::cerr << "[ILP-Gurobi] Error code = " << e.getErrorCode() << ": " << e.getMessage() << std::endl;
        return false;
    } catch (...) {
        if (verbose)
            std::cerr << "[ILP-Gurobi] Unknown exception" << std::endl;
        return false;
    }
}

void recover_s1_skip_comp(NTL::Vec<NTL::ZZ_p> *solution,
                         const unsigned int target_component,
                         const uint8_t *sk,
                         const uint8_t *pk,
                         size_t mlen)
{
    using namespace NTL;

    // storage for up to N equations
    Mat<ZZ_p> eq_sys;
    eq_sys.SetDims(N, N);
    Vec<ZZ_p> eq_sys_line;
    Vec<ZZ_p> eq_sys_b;
    eq_sys_b.SetLength(N);

    uint8_t m[mlen];
    uint8_t sm[mlen + CRYPTO_BYTES];
    size_t smlen;
    unsigned int target_comp_actual = target_component;

    int collected = 0;               // the number of collected equations

    while (collected < N) 
    {
        randombytes(m, mlen);
        crypto_sign_skip(sm, &smlen, m, mlen, sk, target_comp_actual, collected % N);

        faults++;

        try {
            ZZ_p eq_b_comp;
            eq_from_sig(&eq_sys_line, &eq_b_comp, &target_comp_actual,
                        sm, m, target_component, collected % N, pk, mlen);
            eq_sys[collected] = eq_sys_line;
            eq_sys_b[collected] = eq_b_comp;
            collected++;
        } catch (const std::exception &e) {
            continue; // failed to extract equation, try again
        }

        // Try ILP once enough equations have been collected.
        if (collected >= MIN_ROWS_TO_TRY_ILP) {
            std::vector<bool> exact_rows(collected, true);
            NTL::Vec<NTL::ZZ_p> ilp_sol;
            
            bool ok_ilp = false;
            try {
                ok_ilp = solve_ilp_gurobi(ilp_sol,
                                        eq_sys,
                                        eq_sys_b,
                                        collected,
                                        N,
                                        ETA,
                                        exact_rows,
                                        ILP_TIME_LIMIT_MS);
            } catch (...) {
                ok_ilp = false;
            }

            if (ok_ilp) {
                solution->SetLength(N);
                for (int j = 0; j < N; ++j) {
                    (*solution)[j] = ilp_sol[j];
                }
                return;
            }
        }
    }

    throw std::runtime_error("Failed to recover s1 component after collecting N equations and ILP attempts");
}

Vec<ZZ_pX> recover_s1_skip(const uint8_t *sk, const uint8_t *pk, size_t mlen)
{
    Vec<ZZ_p> solution;
    Vec<ZZ_pX> s1_recovered;
    solution.SetLength(N);
    s1_recovered.SetLength(L);

    for (auto j = 0; j < L; j++)
    {
        int safety_stop = 0;
        try_again_1:
        try
        {
            recover_s1_skip_comp(&solution, j, sk, pk, mlen);
        }
        catch(const std::exception& e)
        {
            cout << e.what() << endl;
            safety_stop++;
            if (safety_stop > 10)
            {
                throw logic_error("Something went wrong. s1 component recovery failed 10 times (recover_s1_skip).");
            }
            goto try_again_1;
        }

        for (auto i = 0; i < N; i++)
        {
            SetCoeff(s1_recovered[j], i, solution[i]);
        }
    }

    return s1_recovered;
}

Vec<ZZ_pX> recover_s1_skip_shuff(const uint8_t *sk, const uint8_t *pk, size_t mlen)
{
    Vec<ZZ_pX> s1_recovered;
    s1_recovered.SetLength(L);

    Vec<Mat<ZZ_p>> eq_sys_s;
    eq_sys_s.SetLength(L);
    Vec<Vec<ZZ_p>> eq_sys_b_s;
    eq_sys_b_s.SetLength(L);
    Vec<unsigned int> eq_counts;
    eq_counts.SetLength(L);

    for (int j = 0; j < L; j++)
    {
        eq_sys_s[j].SetDims(N, N);
        eq_sys_b_s[j].SetLength(N);
        eq_counts[j] = 0;
    }

    Vec<ZZ_p> eq_sys_line;
    ZZ_p eq_sys_b_comp;

    uint8_t m[mlen];
    uint8_t sm[mlen + CRYPTO_BYTES];
    size_t smlen;
    unsigned int target_comp_actual = 0;

    std::vector<bool> comp_solved((size_t)L, false);
    int solved_components = 0;

    // Collect equations and try ILP as soon as a component has enough rows.
    for (int i = 0; i < 2 * L * N && solved_components < L; i++)
    {
        randombytes(m, mlen);
        crypto_sign_skip(sm, &smlen, m, mlen, sk, 1, 1);
        faults++;

        try
        {
            eq_from_sig(&eq_sys_line, &eq_sys_b_comp, &target_comp_actual, sm, m, 1, 1, pk, mlen);
        }
        catch (const std::exception&)
        {
            continue;
        }

        if (eq_sys_b_comp == 0)
        {
            continue;
        }
        if (target_comp_actual >= (unsigned int)L)
        {
            continue;
        }
        if (comp_solved[target_comp_actual])
        {
            continue;
        }
        if (eq_counts[target_comp_actual] >= (unsigned int)N)
        {
            continue;
        }

        const unsigned int row = eq_counts[target_comp_actual];
        eq_sys_s[target_comp_actual][row] = eq_sys_line;
        eq_sys_b_s[target_comp_actual][row] = eq_sys_b_comp;
        eq_counts[target_comp_actual] = row + 1;

        const int m_eq = (int)eq_counts[target_comp_actual];
        if (m_eq < MIN_ROWS_TO_TRY_ILP)
        {
            continue;
        }

        std::vector<bool> exact_rows((size_t)m_eq, true);
        Vec<ZZ_p> solution;
        bool ok = false;
        try
        {
            ok = solve_ilp_gurobi(solution,
                                  eq_sys_s[target_comp_actual],
                                  eq_sys_b_s[target_comp_actual],
                                  m_eq,
                                  N,
                                  ETA,
                                  exact_rows,
                                  ILP_TIME_LIMIT_MS);
        }
        catch (...)
        {
            ok = false;
        }

        if (!ok)
        {
            continue;
        }

        for (int j = 0; j < N; j++)
        {
            SetCoeff(s1_recovered[target_comp_actual], j, solution[j]);
        }
        comp_solved[target_comp_actual] = true;
        solved_components++;
    }

    if (solved_components < L)
    {
        throw std::runtime_error("recover_s1_skip_shuff: not all components solved within the iteration bound");
    }

    return s1_recovered;
}