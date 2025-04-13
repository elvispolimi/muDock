from pathlib import Path
import sys

if len(sys.argv) != 3:
    print("Usage: python check_non_bonds.py <autodock_non_bond> <mudock_non_bonf>")
    sys.exit(1)

AUTODOCK_NONBOND = Path(sys.argv[1])
MUDOCK_NONBOND = Path(sys.argv[2])

total = 0
error = 0
with open(MUDOCK_NONBOND, "r") as mudock_non_bond:
    with open(AUTODOCK_NONBOND, "r") as autodock_non_bond:
        mudock_lines = mudock_non_bond.readlines()
        autodock_lines = autodock_non_bond.readlines()
        expected_error = abs(len(mudock_lines) - len(autodock_lines))
        mudock_iter = iter(mudock_lines)
        mudock_line = next(mudock_iter)
        mudock_tokens = mudock_line.split()
        mudock_a1 = int(mudock_tokens[1])
        mudock_a2 = int(mudock_tokens[2])
        mudock_elec = float(mudock_tokens[4])
        mudock_vdw = float(mudock_tokens[5])
        mudock_desolv = float(mudock_tokens[6])
        for autodock_line in autodock_lines:
            tokens = autodock_line.split()
            contrib = float(tokens[3])
            atom_ids = tokens[1].split("-")
            a1 = int(atom_ids[0]) - 1
            a2 = int(atom_ids[1]) - 1
            elec = float(tokens[4])
            vdw = float(tokens[5])
            hb = float(tokens[6])
            desolv = float(tokens[7])
            if a1 != mudock_a1 or a2 != mudock_a2:
                total += contrib
                error += 1
            else:
                if (
                    abs(mudock_elec - elec) > 0.01
                    or abs(mudock_vdw - vdw - hb) > 0.01
                    or abs(mudock_desolv - desolv) > 0.01
                ):
                    raise RuntimeError("Error in intra energies")
                try:
                    mudock_line = next(mudock_iter)
                    mudock_tokens = mudock_line.split()
                    mudock_a1 = int(mudock_tokens[1])
                    mudock_a2 = int(mudock_tokens[2])
                    mudock_elec = float(mudock_tokens[4])
                    mudock_vdw = float(mudock_tokens[5])
                    mudock_desolv = float(mudock_tokens[6])
                except StopIteration:
                    mudock_a1 = -1
                    mudock_a2 = -1
if expected_error != error:
    raise RuntimeError(f"Expected {expected_error} error, but found {error} error")
print(f"total {total}")
print(f"error {error}")
