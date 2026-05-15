import csv
from collections import namedtuple
import sys
import matplotlib.pyplot as plt


mode = sys.argv[1].lower()
OverallTime = {'dup_time': 0, 'inacc_time':0, 'poll_time':0, 'comp_time':0}

filename = "results_combined"+str(mode)+".csv"

PRINT_SUMMARY_TABLE = False
PRINT_TIMES = False
size_to_state={16140044: 'California', 7963659: 'Florida', 5007954: 'Georgia', 1679423: 'Connecticut',1000000: '1 mil', 456000: 'Congress'}
# size_to_state={1000000: '1 mil', 456000: 'Congress', 1679423: 'Connecticut'}
time_results = {}
best_params = {}
best_params_inacc = {}
with open("all_time.csv", "w", newline='') as dupfile:
    writer = csv.writer(dupfile)
    writer.writerow(["size", "margin", "dup time", "inacc time", "poll time", "comp time"])

    mu_list = [.3, .2, .1, .03, .025, .02, .015, .014, .013, .012, .011, .01, .009,.008, .007,.006, .005]
#    mu_list = [0.03, 0.025, 0.02, .019, .018, .017, .016, .015, .014, .013, .012, .011, 0.01, .0075]
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        while True:
            rows = [next(reader, None) for _ in range(3)]
            if None in rows:
                break
            
            try:
                    
                size, margin, k1_max, k2_max, k3_max, in_acc = [float(rows[0][0]), float(rows[0][1]), float(rows[0][4]), float(rows[0][5]), float(rows[0][6]), bool(str(rows[0][14]).lower()=='true')]
                _, _, k1_one, k2_one, k3_one = [float(rows[1][0]), float(rows[1][1]), float(rows[1][4]), float(rows[1][5]), float(rows[1][6])]
                _, _, k1_ten, k2_ten, k3_ten = [float(rows[2][0]), float(rows[2][1]), float(rows[2][4]), float(rows[2][5]), float(rows[2][6])]
            except Exception as e:
                print("Error processing rows:", rows)
                print(e)
                exit(1)
            
            time_comp_max = (35*max(k1_max, k2_max) + 10*k1_max +25*k2_max)/3600
            time_comp_one = (35*max(k1_one, k2_one) + 10*k1_one +25*k2_one)/3600
            time_comp_ten = (35*max(k1_ten, k2_ten) + 10*k1_ten +25*k2_ten)/3600
            time_max = time_comp_max+ k3_max/1538
            time_one = time_comp_one + k3_one/1538
            time_ten = time_comp_ten + k3_ten/1538
            time = min(time_max, time_one, time_ten)
            time_comp = min(time_comp_max, time_comp_one, time_comp_ten)
            if margin not in mu_list:
                continue
            if in_acc:
                OverallTime['inacc_time'] = round(time)
                
                # OverallTime['poll_time'] = round(polling_samples[margin]/60+ size/1538)
                # OverallTime['comp_time'] = round(comp_samples[margin]*(45+25)/3600+size/1538)
                # if size%1000 != 0:
                time_results[(size, margin)] = OverallTime['dup_time']/OverallTime['inacc_time']
                best_params_inacc[(size, margin)] = (k1_max, k2_max, k3_max) if time_max <= time_one and time_max <= time_ten else (k1_one, k2_one, k3_one) if time_one <= time_ten else (k1_ten, k2_ten, k3_ten)
                if PRINT_TIMES:
                    print(size, margin, OverallTime['dup_time'], OverallTime['inacc_time'])
                if time_results[(size, margin)] < 1:
                    print(f"Warning: Duplication time is less than inaccuracy time for size {size} and margin {margin}. Ratio: {time_results[(size, margin)]}")
                # writer.writerow([size, margin, OverallTime['dup_time'], OverallTime['inacc_time'], OverallTime['poll_time'], OverallTime['comp_time']])
            else:
                OverallTime['dup_time'] = round(time)
                best_params[(size, margin)] = (k1_max, k2_max, k3_max) if time_max <= time_one and time_max <= time_ten else (k1_one, k2_one, k3_one) if time_one <= time_ten else (k1_ten, k2_ten, k3_ten)
            
        
        csvfile.close()
    dupfile.close()        


    # Prepare data for plotting
    sizes = sorted(set(size for (size, margin) in time_results.keys()))
    margins = sorted(set(margin for (_, margin) in time_results.keys()))

    plt.figure(figsize=(10, 6))

    for size in sizes:

        y = [time_results.get((size, margin), None) for margin in margins]
        plt.plot(margins[:len(y)], y, marker='o', label=f'Size: {size_to_state.get(size, size)}')

    plt.xlabel("Margin")
    plt.ylabel("Time Ratio")
    if mode in ["comparison", 'polling']:
        plt.yscale("log")
    plt.title("Time Results by Size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("time_results_by_size"+str(mode)+".png")
    # plt.show()
    
    # Save data as CSV
    with open("time_results_by_size"+str(mode)+".csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Margin"] + [size_to_state.get(size, size) for size in sizes])
        for margin in margins:
            row = [margin] + [time_results.get((size, margin), None) for size in sizes]
            writer.writerow(row)

margins = [.3, .2, .1, .03, .02, .01, .005]


size_lookup = {456000: 'Cong.', 1000000: '1 M', 1679423: 'CT', 5007954: 'GA', 7963659: 'FL', 16140044: 'CA'}

if PRINT_SUMMARY_TABLE:
    # Generate LaTeX table
    print("\\begin{table*}[t]")
    print("\\centering")
    print("\\begin{tabular}{|cc|ccc|ccc|}")
    print("\\hline")
    print("& & & & \\multicolumn{3}{c|}{Accurate} & \\multicolumn{3}{c}{Inaccurate} \\\\")
    print("$N$ & $\\mu$ & $k_{\\mathrm{dup}}$ & $k_{\\mathrm{comp}}$ & $k_{\\mathrm{size}}$ &")
    print("$k_{\\mathrm{dup}}$ & $k_{\\mathrm{comp}}$ & $k_{\\mathrm{size}}$ \\\\")
    print("\\hline")

    for size in sizes:
        for margin in sorted(margins):
            if (size, margin) in sorted(best_params) and (size, margin) in sorted(best_params_inacc):
                k1_acc, k2_acc, k3_acc = best_params[(size, margin)]
                k1_inacc, k2_inacc, k3_inacc = best_params_inacc[(size, margin)]
                if margin ==.005:
                    print("\\multirow{7}{*}{" + size_lookup.get(size) + "} ", end="")
                if mode != "direct":
                    if mode == "polling":
                        print(f"& {margin*100:.2f}\\% & {int(k2_acc/1000):d} & {int(size/1000):d} & {int(k2_inacc/1000):d} & {int(k3_inacc/1000):d} \\\\")
                    else:
                        print(f"& {margin*100:.2f}\\% & {int(k2_acc/1):d} & {int(size/1000):d} & {int(k2_inacc/1):d} & {int(k3_inacc/1000):d} \\\\")
                else:
                    print(f"& {margin*100:.2f}\\% & {int(k1_acc/1000):d} & {int(k2_acc/1000):d} & {int(size/1000):d} & {int(k1_inacc/1000):d} & {int(k2_inacc/1000):d} & {int(k3_inacc/1000):d} \\\\")
                if margin == .3:
                    print("\\hline")

    print("\\end{tabular}")
    print("\\end{table*}")