import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f:
    lines = f.readlines()

with open(output_file, "w") as out:
    seq = ""
    header = ""

    for line in lines:
        line = line.strip()

        if line.startswith(">P1;"):
            if seq:
                out.write(f">{header}\n{seq}\n")
                seq = ""
            header = line.replace(">P1;", "").strip()

        elif line.startswith("structure") or line.startswith("sequence"):
            continue

        elif not line:
            continue

        else:
            clean_line = line.replace("*", "").replace(" ", "").upper()
            seq += clean_line

    if seq:
        out.write(f">{header}\n{seq}\n")
