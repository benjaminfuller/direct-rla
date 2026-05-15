from math import log, ceil, sqrt, exp
import numpy as np
import sys
from scipy.stats import hypergeom
from scipy.stats import nchypergeom_fisher as fisher
import csv
from functools import lru_cache
from Minerva_Sample_size import round_size_approx

# ===================== CONSTANTS =====================
average_batch_size = 900
inaccConstant = .10
inaccConstantSmall = .001
gamma = 1.1
num_trials = 1000

# ===================== LOAD CSV =====================
size_per_risk = {}
all_alpha = set()
all_mu = set()
size_items = []
o1rate = .001
u1rate = .001
o2rate = .0001
u2rate = .0001
gamma = 1.1

with open('ComparisonSizes.csv', 'r', newline='') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        alpha = float(row[1])
        mu = float(row[0]) / 100
        size = int(row[2])

        all_alpha.add(alpha)
        all_mu.add(mu)

        if (alpha, mu) not in size_per_risk or size < size_per_risk[(alpha, mu)]:
            size_per_risk[(alpha, mu)] = size

    print(f"Loaded {len(size_per_risk)} entries from CSV. Unique alpha: {len(all_alpha)}, Unique mu: {len(all_mu)}")

size_items = [(a, m, s) for (a, m), s in size_per_risk.items()]

def compute_kaplan_markov_size(alpha, mu):
    outcomes = np.array([-2, -1, 0, 1, 2])
    probabilities = np.array([u2rate, u1rate, 1 - o1rate - o2rate - u1rate - u2rate, o1rate, o2rate])
    base_factor = (1 - (mu / (2 * gamma)))
    k = 0
    observedrisk = 1.0
    while k <= 500000:
        samples = np.random.choice(outcomes, size=10000, p=probabilities)
        factors = base_factor / (1 - (samples / (2 * gamma)))
        cumprod = np.cumprod(factors)
        below_alpha = np.where(observedrisk * cumprod < alpha)[0]
        if below_alpha.size > 0:
            return k + below_alpha[0] + 1
        observedrisk *= cumprod[-1]
        k += len(samples)
    return 2**32

def run_sim_missing(alpha, mu):
    k_array = np.array([compute_kaplan_markov_size(alpha, mu) for _ in range(num_trials)])
    return np.percentile(k_array, 90)

def lookup_comparison_size(alpha_2, mu_2):
    # print(mu_2, alpha_2, flush=True)
    if mu_2 <=0:
        return None
    best = None
    for a, m, s in size_items:
        if a <= alpha_2 and m <= mu_2 and (m >  .02 or (abs(a-alpha_2)<= .002 and abs(m-mu_2) <= .001)):
        # if a <= alpha_2 and m <= mu_2:
            if best is None or s < best:
                best = s
    if best is None:
        k_2 = run_sim_missing(alpha_2, mu_2)    
        size_items.append((alpha_2, mu_2, k_2))
        size_per_risk[(alpha_2, mu_2)] = k_2
        best = k_2

        if len(size_items) %100 == 0:
            with open('ComparisonSizes.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                for alpha_2, mu_2 in size_per_risk:
                    k_2 = size_per_risk[(alpha_2, mu_2)]
                    writer.writerow([round(mu_2 * 100, 5), round(alpha_2, 4), int(k_2)])
            print(f"Size items length: {len(size_items)}", flush=True)
    return best

# ===================== PRECOMPUTATION =====================
def compute_duplicate_risk(N, d, t):
    temp0 = -t/N
    temp = exp(temp0)
    temp2 = -d*(1-temp)**2
    return 2*exp(temp2) 

# ===================== SEARCH =====================

def find_k_from_risk(N, N_eff, duplicates, risk_array, alpha, k_candidates):
    lo, hi = 0, len(k_candidates) - 1
    best_k = None
    max_k =  int(min(1000000, N)) 
    while lo <= hi:
        mid = (lo + hi) // 2
        k = int(k_candidates[mid])
        if (N_eff, duplicates, k)  not in risk_array:
            risk_array[N_eff, duplicates, k] = compute_duplicate_risk(N_eff, duplicates, k)
        if risk_array[N_eff, duplicates, k] <= alpha:
            best_k = k
            hi = mid - 1
        else:
            lo = mid + 1
    return best_k

# ===================== MAIN =====================

mu_list = [
    .30, .20, .10, 0.03, 0.025, 0.02, .019, .018, .017, .016, .015,
    .014, .013, .012, .011, 0.01, .009, .008,
    .0075, .007, .006, 0.005, 0.004, 0.003, 0.002
]

def run_simulation(output_csv, N, mode):

    N_eff = int(N * (1 + inaccConstant))
    max_k = int(min(1000000, N))
    k_candidates = np.concatenate([np.arange(100, 1000, 50), np.arange(1000, max_k, 500)])

    if mode in  ['direct']:
        duplicates_values = list()
        acc_risks = {}
        inacc_risks = {}


    with open(output_csv, 'w', newline='') as out_file:
        writer = csv.writer(out_file)

        writer.writerow([
            "N", "mu", "alpha", "method",
            "k_1", "k_2", "k_3",
            "alpha1","alpha2","alpha3",
            "mu1","mu2","mu3",
            "mode", "INACCURATE_MANIFEST"
        ])

        for mu_original in mu_list:
            for INACCURATE_MANIFEST in [False, True]:

                if INACCURATE_MANIFEST:
                    size_array = np.arange(.95, .01, -.005)
                else:
                    size_array = [0]

                best = {}

                if mode in ['comparison', 'polling']:
                    alpha_var_range = [0]
                    mu_var_range = [0]
                else:
                    alpha_var_range = np.arange(0.05, 1, 0.025)
                    mu_var_range = np.arange(0.05, .6, 0.025)
                
                for alpha_var in alpha_var_range:
                    for mu_var in mu_var_range:

                        for portion_to_size in size_array:

                            mu = mu_original
                            alpha = 0.05
                            k_3 = N

                            if INACCURATE_MANIFEST:
                                p = portion_to_size

                                if p * mu <= inaccConstantSmall:
                                    p = 1.5 * inaccConstantSmall / mu
                                    if p > .95:
                                        continue

                                alpha_3 = p * alpha
                                mu_3 = p * mu
                                alpha = (1 - p) * alpha

                                p_size = (mu_3 - inaccConstantSmall) / (inaccConstant - inaccConstantSmall)

                                tau = 1 / (
                                    (1 - p_size) / (1 + inaccConstantSmall)
                                    + p_size / (1 + inaccConstant)
                                ) - 1

                                total_variation = (
                                    inaccConstantSmall / (2 + inaccConstantSmall)
                                    + (1 + inaccConstantSmall) * p_size * inaccConstant
                                      / (1 + inaccConstant * p_size)
                                )

                                epsilon_p = (1-p_size)*inaccConstantSmall + p_size*inaccConstant
                                N_eff = int(N*(1+epsilon_p))
                                eta = (1+inaccConstant)*(1+epsilon_p)-1
                                if eta > sqrt(2)-1:
                                    continue

                                if mode == 'polling':
                                    mu = mu_original - tau - 2 * total_variation
                                    mu_3 = mu_original - tau - 2 * total_variation
                                elif mode == 'direct':
                                    mu = mu_original - tau - 4 * total_variation
                                    mu_3 = mu_original - tau - 4 * total_variation
                                else:
                                    mu = mu_original - epsilon_p
                                    mu_3 = epsilon_p
                                k_3 = ceil(-(1 + inaccConstantSmall) * log(alpha_3) / p_size)

                                if k_3 < 0:
                                    continue

                                k_3 = average_batch_size * ceil(k_3)

                                if k_3 >= 3 * N / 4 or mu <= 0:
                                    k_3 = N
                                    alpha_3 = 0
                                    mu_3 = 0
                                    alpha = 0.05
                                    mu = mu_original
                                    tau = 0
                            else:
                                alpha_3 = 0
                                mu_3 = 0
                                tau = 0

                            # ---- K1 ----
                            mu_1 = mu * mu_var
                            alpha_1 = alpha * alpha_var

                            if mode == 'direct':
                                duplicates = int(mu_1 * N_eff)
                                if INACCURATE_MANIFEST and k_3 < N:
                                    k_1 = find_k_from_risk(N, N_eff, duplicates,inacc_risks,alpha_1,k_candidates)
                                else:
                                    k_1 = find_k_from_risk(N, N, duplicates,acc_risks,alpha_1,k_candidates)
            

                                if k_1 is None:
                                    continue
                            else:
                                k_1 = 0
                                alpha_1 =0
                                
                            alpha_2 = round(alpha * (1 - alpha_var), 7)
                            mu_2 = round(mu * (1 - 2*mu_var), 7)

                            if INACCURATE_MANIFEST and mode in ['polling', 'direct']:
                                mu_2 *= (1 + tau)
                                mu_2 = round(mu_2, 7)

                            if mode == 'polling':
                                k_2 = round_size_approx(mu_2, alpha_2, 0.9)
                            else:
                                k_2 = lookup_comparison_size(alpha_2, mu_2)
                            if k_2 is None:
                                continue
                            
                            candidates = {
                                'max': max(k_1, k_2),
                                'one': 10*k_1 + 25*k_2 + 25*max(k_1, k_2)+ k_3/1538*3600,
                                'ten': 10*k_1 + 62*k_2 + 25*max(k_1, k_2)+ k_3/1538*3600
                            }
                            for key, val in candidates.items():
                                if key not in best or val < best[key][0]:
                                    best[key] = (val, k_1, k_2, k_3, mu_1, mu_2, mu_3, alpha_1, alpha_2, alpha_3)
                                    # print(best[key][0], val, k_1, k_2, k_3, mu_1, mu_2, mu_3, alpha_1, alpha_2, alpha_3)

                for key in ['max', 'one', 'ten']:
                    if key in best:
                        _, k1, k2, k3, mu1, mu2, mu3, a1, a2, a3 = best[key]

                        writer.writerow([
                            N, round(mu_original, 4), 0.05, key,
                            int(k1), int(k2), int(min(N, k3)),
                            round(a1, 4), round(a2, 4), round(a3, 4),
                            round(mu1, 7), round(mu2, 7), round(mu3, 7),
                            mode, INACCURATE_MANIFEST
                        ])
                        # print(N, round(mu_original, 4), key, int(k1), int(k2), int(k3), round(a1, 4), round(a2, 4), round(a3, 4), INACCURATE_MANIFEST,flush=True)
                        out_file.flush()
    


# ===================== RUN =====================

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    mode = sys.argv[2].lower() if len(sys.argv) > 2 else 'direct'

    if mode not in ['direct', 'polling', 'comparison']:
        raise ValueError("Mode must be one of: direct, polling, comparison ", mode)

    run_simulation(f'results_combined{N}_{mode}.csv', N, mode)