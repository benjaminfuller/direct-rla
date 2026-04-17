import csv

# These were directly pulled from Minerva approximation but can be replaced as needed
MINERVA = {0.03: 9844, 0.025: 14176, 0.02: 22152, 0.015: 39384, 0.01: 88619, 0.0075: 157547, 0.005: 354486}

mu_list = [.03, .025, .02, .015, .01, .0075, .005]
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
                k_dup = int(r["k_1"]) *900/int(r["N"])
                k_sam = int(r["k_2"]) *900/int(r["N"])

                a_dup = int(float(r["alpha1"]) /.05 * 100)
                a_sam = int(float(r["alpha2"]) /.05 * 100)

                minerva = MINERVA.get(mu, "") *900/int(r["N"])
                rows.append((size, mu, k_dup, a_dup, k_sam, a_sam, minerva))

    rows.sort(key=lambda x: (x[0], x[1]))

    for N, mu, k1, a1, k2, a2, m in rows:
        print(
            f"& {format_margin(mu):.2f} & "
            f"{k1:.3f} & {a1} & "
            f"{k2:.3f} & {a2} & {int(round(m,0))} \\\\"
        )

if __name__ == "__main__":
    main()
