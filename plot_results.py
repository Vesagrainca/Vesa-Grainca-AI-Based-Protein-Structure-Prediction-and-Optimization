import csv
from pathlib import Path
from statistics import mean
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
with (ROOT / 'experiment_results.csv').open(newline='', encoding='utf-8-sig') as file:
    rows = list(csv.DictReader(file))

output = ROOT / 'protein_graphs'
output.mkdir(exist_ok=True)
lengths = sorted({int(r['length']) for r in rows})
algorithms = ['Genetic Algorithm', 'Simulated Annealing', 'Hybrid Algorithm']
metrics = [('energy', 'Mean best energy', 'Energy comparison'), ('time', 'Mean execution time (s)', 'Execution time comparison'), ('memory', 'Mean peak Python memory (KiB)', 'Memory comparison'), ('valid', 'Mean valid evaluations (%)', 'Valid evaluation rate'), ('pruning', 'Mean pruned moves (%)', 'Pruning comparison')]
for metric, ylabel, title in metrics:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for algorithm in algorithms:
        values = []
        for n in lengths:
            group = [r for r in rows if int(r['length']) == n and r['algorithm'] == algorithm]
            values.append(mean(float(r[metric]) / float(r['explored']) * 100 if metric == 'valid' else float(r[metric]) for r in group))
        ax.plot(lengths, values, marker='o', label=algorithm)
    ax.set(xlabel='Protein length', ylabel=ylabel, title=title, xticks=lengths)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / f'{metric}_comparison.png', dpi=180)
    plt.close(fig)
print(f'Graphs saved in {output}')
