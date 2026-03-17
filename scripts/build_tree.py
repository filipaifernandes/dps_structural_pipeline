import sys
import subprocess
import os

alignment_file = sys.argv[1]
output_tree = sys.argv[2]

os.makedirs(os.path.dirname(output_tree), exist_ok=True)

cmd = f"FastTree {alignment_file} > {output_tree}"
subprocess.run(cmd, shell=True, check=True)

print(f"Tree written to {output_tree}")
