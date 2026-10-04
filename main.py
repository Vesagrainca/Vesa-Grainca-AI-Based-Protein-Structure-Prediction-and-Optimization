class Protein:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def __str__(self):
        return f"Protein(sequence='{self.sequence}', length={self.length})"


def directions_to_coordinates(directions):
    coordinates = [(0, 0)]

    x, y = 0, 0

    for direction in directions:
        if direction == "R":
            x += 1
        elif direction == "L":
            x -= 1
        elif direction == "U":
            y += 1
        elif direction == "D":
            y -= 1
        else:
            raise ValueError(f"Invalid direction: {direction}")

        coordinates.append((x, y))

    return coordinates


def manhattan_distance(p1, p2):
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])


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
        return -2

    if {type1, type2} == {"+", "-"}:
        return -2

    if type1 == "P" and type2 == "P":
        return -1

    if type1 == type2 and type1 in {"+", "-"}:
        return 1

    return 0


def calculate_energy(protein, coordinates):
    if not is_valid_structure(coordinates):
        return float("inf")

    energy = 0

    for i in range(protein.length):
        for j in range(i + 2, protein.length):
            if manhattan_distance(coordinates[i], coordinates[j]) == 1:
                type1 = amino_acid_type(protein.sequence[i])
                type2 = amino_acid_type(protein.sequence[j])

                energy += interaction_energy(type1, type2)

    return energy


# -------------------------
# TESTING CODE
# -------------------------

protein = Protein("ACDEF")

directions = ["R", "R", "U", "R"]

coordinates = directions_to_coordinates(directions)

print(protein)
print("Directions:", directions)
print("Coordinates:", coordinates)
print("Valid structure:", is_valid_structure(coordinates))

energy = calculate_energy(protein, coordinates)

print("Energy:", energy)                                             