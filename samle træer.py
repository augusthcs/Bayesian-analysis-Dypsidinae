import os

# Directory containing Newick files
directory = "C:\\Users\\augus\\SKOLe\\5. semester\\Bachelor\\contrees"

# File to store combined trees
output_file = "input_trees.newick"

# Open output file for writing
with open(output_file, "w") as outfile:
    for filename in os.listdir(directory):
        if filename.endswith(".contree"):  # Check for Newick files
            with open(os.path.join(directory, filename), "r") as infile:
                tree = infile.read().strip()
                outfile.write(tree + "\n")
