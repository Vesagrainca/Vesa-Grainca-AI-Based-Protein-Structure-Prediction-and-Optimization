import random
import math
import time
import tracemalloc
import csv
import matplotlib.pyplot as plt

statistics = {
    "explored": 0,
    "valid": 0,
    "invalid": 0,
    "checked_moves": 0,
    "pruned_moves": 0,
    "best_energy": float("inf"),
    "convergence": []
}
class Protein:
    def __init__(self, sequence):
        sequence = sequence.strip().upper()

        if not sequence:
            raise ValueError("Protein sequence cannot be empty")

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

        if any(aa not in valid_amino_acids for aa in sequence):
            raise ValueError("Invalid amino acid in sequence")

        self.sequence = sequence
        self.length = len(sequence)

    def __str__(self):
        return f"Protein(sequence='{self.sequence}', length={self.length})"
def directions_to_coordinates(directions):
    coordinates = [(0, 0, 0)]

    x, y, z = 0, 0, 0

    for direction in directions:
        if direction == "R":
            x += 1
        elif direction == "L":
            x -= 1
        elif direction == "U":
            y += 1
        elif direction == "D":
            y -= 1
        elif direction == "F":
            z += 1
        elif direction == "B":
            z -= 1
        else:
            raise ValueError(f"Invalid direction: {direction}")

        coordinates.append((x, y, z))

    return coordinates


def manhattan_distance(p1, p2):
    return (
        abs(p1[0] - p2[0])
        + abs(p1[1] - p2[1])
        + abs(p1[2] - p2[2])
    )


def is_valid_structure(coordinates):
    # Two amino acids cannot occupy the same position
    if len(coordinates) != len(set(coordinates)):
        return False

    # Consecutive amino acids must remain next to each other
    for i in range(len(coordinates) - 1):
        if manhattan_distance(coordinates[i], coordinates[i + 1]) != 1:
            return False

    return True


# Amino-acid categories
ENERGY_PARAMETERS = {
    "hydrophobic": -2,
    "opposite_charge": -2,
    "polar": -1,
    "same_charge": 1
}
HYDROPHOBIC = set("AVILMFWY")
POLAR = set("STNQCGP")
POSITIVE = set("KRH")
NEGATIVE = set("DE")


def amino_acid_type(amino_acid):
    if amino_acid in HYDROPHOBIC:
        return "H"

    if amino_acid in POLAR:
        return "P"

    if amino_acid in POSITIVE:
        return "+"

    if amino_acid in NEGATIVE:
        return "-"

    return "P"


def interaction_energy(type1, type2):
    if type1 == "H" and type2 == "H":
        return ENERGY_PARAMETERS["hydrophobic"]

    if {type1, type2} == {"+", "-"}:
        return ENERGY_PARAMETERS["opposite_charge"]

    if type1 == "P" and type2 == "P":
        return ENERGY_PARAMETERS["polar"]

    if type1 == type2 and type1 in {"+", "-"}:
        return ENERGY_PARAMETERS["same_charge"]

    return 0

def calculate_energy(protein, coordinates):
    if len(coordinates) != protein.length:
        return float("inf")

    if not is_valid_structure(coordinates):
        return float("inf")

    energy = 0

    positions = {
        position: i
        for i, position in enumerate(coordinates)
    }

    neighbors = [
        (1, 0, 0), (-1, 0, 0),
        (0, 1, 0), (0, -1, 0),
        (0, 0, 1), (0, 0, -1)
    ]

    for i, (x, y, z) in enumerate(coordinates):
        for dx, dy, dz in neighbors:

            neighbor = (x + dx, y + dy, z + dz)
            j = positions.get(neighbor)

            if j is not None and j > i + 1:
                type1 = amino_acid_type(protein.sequence[i])
                type2 = amino_acid_type(protein.sequence[j])

                energy += interaction_energy(type1, type2)

    return energy
    return energy
def generate_random_structure(protein):
    possible_directions = {
        "R": (1, 0, 0),
        "L": (-1, 0, 0),
        "U": (0, 1, 0),
        "D": (0, -1, 0),
        "F": (0, 0, 1),
        "B": (0, 0, -1)
    }

    directions = []
    occupied = {(0, 0, 0)}

    x, y, z = 0, 0, 0

    for i in range(protein.length - 1):
        valid_directions = []

        for direction, (dx, dy, dz) in possible_directions.items():
            if i == 0 and direction != "R":
                continue

            new_position = (x + dx, y + dy, z + dz)

            statistics["checked_moves"] += 1

            if new_position not in occupied:
                valid_directions.append((direction, new_position))
            else:
                statistics["pruned_moves"] += 1

        if not valid_directions:
            return ["R"] * (protein.length - 1)

        direction, new_position = random.choice(valid_directions)

        directions.append(direction)
        occupied.add(new_position)

        x, y, z = new_position

    return directions
def initialize_population(protein, population_size):
    population = []

    for i in range(population_size):
        individual = generate_random_structure(protein)
        population.append(individual)

    return population
def evaluate_individual(protein, individual):
    statistics["explored"] += 1

    coordinates = directions_to_coordinates(individual)
    energy = calculate_energy(protein, coordinates)

    if energy == float("inf"):
        statistics["invalid"] += 1
    else:
        statistics["valid"] += 1

    if energy < statistics["best_energy"]:
        statistics["best_energy"] = energy

    statistics["convergence"].append(
        statistics["best_energy"]
    )

    return energy

    return sorted_population[0], sorted_population[1]
def crossover(parent1, parent2):
    if len(parent1) < 2:
        return parent1.copy()

    crossover_point = random.randint(1, len(parent1) - 1)

    child = parent1[:crossover_point] + parent2[crossover_point:]

    return child
def mutate(individual):
    mutated = individual.copy()

    if len(mutated) <= 1:
        return mutated

    possible_directions = ["R", "L", "U", "D", "F", "B"]

    index = random.randint(1, len(mutated) - 1)

    possible_directions.remove(mutated[index])

    mutated[index] = random.choice(possible_directions)

    return mutated
def genetic_algorithm(protein, population_size=50,
                      budget=10000, mutation_rate=0.2):

    if population_size < 2 or budget < population_size:
        raise ValueError("Invalid population size or budget")

    start_evaluations = statistics["explored"]

    population = initialize_population(protein, population_size)

    scored_population = [
        (evaluate_individual(protein, individual), individual)
        for individual in population
    ]

    best_energy, best_individual = min(
        scored_population, key=lambda item: item[0]
    )

    while statistics["explored"] - start_evaluations < budget:

        scored_population.sort(key=lambda item: item[0])

        parents = scored_population[:max(2, population_size // 2)]
        new_population = [scored_population[0]]

        while len(new_population) < population_size:

            if statistics["explored"] - start_evaluations >= budget:
                break

            parent1, parent2 = random.sample(parents, 2)

            child = crossover(parent1[1], parent2[1])

            if random.random() < mutation_rate:
                child = mutate(child)

            energy = evaluate_individual(protein, child)

            if energy == float("inf"):
                if statistics["explored"] - start_evaluations >= budget:
                    break

                child = generate_random_structure(protein)
                energy = evaluate_individual(protein, child)

            new_population.append((energy, child))

            if energy < best_energy:
                best_energy = energy
                best_individual = child.copy()

        scored_population = new_population

    return best_individual, best_energy
def simulated_annealing(protein, budget=10000,
                        temperature=100, cooling_rate=0.995,
                        initial_structure=None):

    if budget < 1:
        raise ValueError("Budget must be positive")

    if initial_structure is None:
        current = generate_random_structure(protein)
    else:
        current = initial_structure.copy()

    current_energy = evaluate_individual(protein, current)

    best = current.copy()
    best_energy = current_energy

    for i in range(budget - 1):

        neighbor = mutate(current)
        neighbor_energy = evaluate_individual(protein, neighbor)

        if neighbor_energy != float("inf"):

            delta = neighbor_energy - current_energy

            if delta <= 0 or random.random() < math.exp(
                -delta / temperature
            ):
                current = neighbor
                current_energy = neighbor_energy

                if current_energy < best_energy:
                    best = current.copy()
                    best_energy = current_energy

        temperature = max(temperature * cooling_rate, 0.000001)

    return best, best_energy
def hybrid_algorithm(protein, budget=10000):

    ga_budget = budget // 2
    sa_budget = budget - ga_budget

    ga_structure, ga_energy = genetic_algorithm(
        protein,
        budget=ga_budget
    )

    best_structure, best_energy = simulated_annealing(
        protein,
        budget=sa_budget,
        initial_structure=ga_structure
    )

    return best_structure, best_energy
def measure_performance(algorithm, protein):

    def reset_statistics():
        for key in statistics:
            statistics[key] = 0

        statistics["best_energy"] = float("inf")
        statistics["convergence"] = []

    reset_statistics()

    random_state = random.getstate()

    start_time = time.perf_counter()

    structure, energy = algorithm(protein)

    execution_time = time.perf_counter() - start_time

    algorithm_stats = statistics.copy()

    random.setstate(random_state)

    reset_statistics()

    tracemalloc.start()

    algorithm(protein)

    current, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()

    memory_usage = peak / 1024

    return structure, energy, execution_time, memory_usage, algorithm_stats
def visualize_protein(protein, directions, title):
    coordinates = directions_to_coordinates(directions)

    x, y, z = zip(*coordinates)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    colors = {
        "H": "orange",
        "P": "green",
        "+": "blue",
        "-": "red"
    }

    ax.plot(x, y, z, color="gray", linewidth=2)

    for i, amino_acid in enumerate(protein.sequence):
        color = colors[amino_acid_type(amino_acid)]

        ax.scatter(x[i], y[i], z[i], color=color, s=150)
        ax.text(x[i], y[i], z[i], amino_acid)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title)
    ax.set_box_aspect((1, 1, 1))

    plt.show()
def plot_convergence(convergence_results):
    plt.figure(figsize=(10, 6))

    for name, history in convergence_results.items():
        plt.plot(
            range(1, len(history) + 1),
            history,
            label=name
        )

    plt.xlabel("Number of evaluations")
    plt.ylabel("Best energy")
    plt.title("Algorithm Convergence Comparison")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.show()
def run_experiments():
    lengths = [10, 20, 30, 50, 100]
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    algorithms = {
        "Genetic Algorithm": genetic_algorithm,
        "Simulated Annealing": simulated_annealing,
        "Hybrid Algorithm": hybrid_algorithm
    }

    results = []

    for length in lengths:
        rng = random.Random(42 + length)
        sequence = "".join(
            rng.choice(amino_acids) for _ in range(length)
        )

        protein = Protein(sequence)

        for trial in range(3):
            for name, algorithm in algorithms.items():

                random.seed(100 + trial)

                structure, energy, elapsed, memory, stats = measure_performance(
                    algorithm, protein
                )

                checked = stats["checked_moves"]
                pruned = stats["pruned_moves"]

                results.append({
                    "length": length,
                    "sequence": sequence,
                    "trial": trial + 1,
                    "algorithm": name,
                    "energy": energy,
                    "time": elapsed,
                    "memory": memory,
                    "explored": stats["explored"],
                    "valid": stats["valid"],
                    "invalid": stats["invalid"],
                    "structural_validity": is_valid_structure(
                        directions_to_coordinates(structure)
                    ),
                    "pruning": pruned / checked * 100 if checked else 0
                })

                print(
                    f"Length: {length} | Trial: {trial + 1} | "
                    f"{name} | Energy: {energy}"
                )

    with open("experiment_results.csv", "w", newline="") as file:
        writer = csv.DictWriter(
            file, fieldnames=results[0].keys()
        )
        writer.writeheader()
        writer.writerows(results)

    print("\nExperiments completed successfully!")

    return results
# -------------------------
# TESTING CODE
# -------------------------

# -------------------------
# TESTING CODE
# -------------------------
def configure_energy():
    choice = input("Customize energy parameters? (y/n): ")

    if choice.lower() == "y":
        new_parameters = {}

        for key, default in ENERGY_PARAMETERS.items():
            value = input(f"{key} [{default}]: ")

            if value.strip():
                new_parameters[key] = int(value)
            else:
                new_parameters[key] = default

        ENERGY_PARAMETERS.update(new_parameters)

def main():
    print("\nAI-Based Protein Structure Prediction")
    print("1. Genetic Algorithm")
    print("2. Simulated Annealing")
    print("3. Hybrid Algorithm")
    print("4. Run experiments")
    print("5. Exit")

    choice = input("Select an option: ")

    algorithms = {
        "1": genetic_algorithm,
        "2": simulated_annealing,
        "3": hybrid_algorithm
    }

    if choice == "4":
        run_experiments()
        return

    if choice == "5":
        return

    if choice not in algorithms:
        print("Invalid option")
        return

    try:
        sequence = input("Enter amino acid sequence: ")
        protein = Protein(sequence)

        budget = int(input("Enter evaluation budget: "))

        if budget < 100:
            raise ValueError("Minimum budget is 100")
        configure_energy()

        algorithm = algorithms[choice]

        structure, energy, elapsed, memory, stats = measure_performance(
            lambda p: algorithm(p, budget=budget), protein
        )

        print("\n--- RESULTS ---")
        print("Protein length:", protein.length)
        print("Best structure:", structure)
        print("Coordinates:", directions_to_coordinates(structure))
        print("Best energy:", energy)
        print("Valid structure:", is_valid_structure(
            directions_to_coordinates(structure)
        ))
        print(f"Execution time: {elapsed:.4f} seconds")
        print(f"Memory usage: {memory:.2f} KB")
        print("Explored:", stats["explored"])
        print("Valid:", stats["valid"])
        print("Invalid:", stats["invalid"])

        checked = stats["checked_moves"]
        pruned = stats["pruned_moves"]
        percentage = pruned / checked * 100 if checked else 0

        print(f"Pruning percentage: {percentage:.2f}%")

        visualize_protein(protein, structure, f"Energy: {energy}")
        plot_convergence({algorithm.__name__: stats["convergence"]})

    except ValueError as error:
        print("Error:", error)


if __name__ == "__main__":
    main()