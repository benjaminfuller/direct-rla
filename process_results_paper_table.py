import csv

# Fixed Minerva values (from your original table)
MINERVA = {0.03: 9844, 0.025: 14176, 0.02: 22152, 0.015: 39384, 0.01: 88619, 
           .0095: 98193, .0090: 109407, 
           .0085: 122657, .0080: 138469, 
           0.0075: 157547, .0070: 180858,
           .0060: 246170, 0.005: 354486,
           .0040: 553886, .0030: 984689, 
           .0020: 2215553}

PER_BATCH = False # Set to True to scale values per batch, False to keep original values
mu_list = [.01, .009, .008, .007, .006, .005, .004, .003, .002]
def format_margin(mu):
    return round(mu * 100, 2)

def main():
    input_file = "results_combineddirect.csv"

    rows = []

    with open(input_file, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r["method"] == "max" and r["INACCURATE_MANIFEST"] == "False":
                mu = float(r["mu"])
                if mu not in mu_list:
                    continue

                # divide by 1000
                size = int(r["N"]) 
                if PER_BATCH:
                    k_dup = int(r["k_1"]) *900/int(r["N"])
                    k_sam = int(r["k_2"]) *900/int(r["N"])
                    minerva = MINERVA.get(mu, "") *900/int(r["N"])

                else:
                    k_dup = int(r["k_1"])/1000
                    k_sam = int(r["k_2"])/1000
                    minerva = MINERVA.get(mu, "")/1000


                a_dup = int(round(float(r["alpha1"]) /.05 * 100))
                a_sam = int(round(float(r["alpha2"]) /.05 * 100))


                rows.append((size, mu, k_dup, a_dup, k_sam, a_sam, minerva))

    rows.sort(key=lambda x: (x[0], x[1]))

    N_temp = None
    for N, mu, k1, a1, k2, a2, m in rows:
        if N != N_temp:
            if N_temp is not None:
                print()  # Newline between different N values
            print(f"N={N} ", end='')
            N_temp = N
        print(f"({format_margin(mu)},{int(k1)}) ", end='')
        # print(
        #     f"& {format_margin(mu):.2f} & "
        #     f"{k1:.0f} & {a1} & "
        #     f"{k2:.0f} & {a2} & {int(round(m,0))} \\\\ \\cline{{2-7}}"
        # )
    

if __name__ == "__main__":
    main()
