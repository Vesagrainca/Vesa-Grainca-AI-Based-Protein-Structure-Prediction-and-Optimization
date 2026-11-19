from main import *

# Test 1: Valid protein
protein = Protein("AAAA")
assert protein.length == 4

# Test 2: Valid structure
coordinates = directions_to_coordinates(["R", "U", "L"])
assert is_valid_structure(coordinates)

# Test 3: Invalid structure
coordinates = directions_to_coordinates(["R", "L"])
assert not is_valid_structure(coordinates)

# Test 4: Energy calculation
protein = Protein("AAAA")
coordinates = directions_to_coordinates(["R", "U", "L"])
assert calculate_energy(protein, coordinates) == -2

# Test 5: Algorithm validation
protein = Protein("ACDEFGHIKLMNPQRSTVWY")

for algorithm in [
    genetic_algorithm,
    simulated_annealing,
    hybrid_algorithm
]:
    structure, energy = algorithm(protein, budget=1000)
    coordinates = directions_to_coordinates(structure)

    assert len(coordinates) == protein.length
    assert is_valid_structure(coordinates)
    assert energy == calculate_energy(protein, coordinates)

print("All tests passed successfully!")