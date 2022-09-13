# Tests results of new hairy graph generation algo for equality

import os


def get_g6_files(path):
    for file in os.listdir(path):
        fullfile = os.path.join(path, file)
        if file.endswith(".g6") and os.path.isfile(fullfile):
            yield file


for path1, path2 in [("gh_data/data/hairy/even_edges_odd_hairs", "gh_data/data/hairy3/even_edges_odd_hairs"),
                     ("gh_data/data/hairy/even_edges_even_hairs",
                         "gh_data/data/hairy3/even_edges_even_hairs"),
                     ("gh_data/data/hairy/odd_edges_odd_hairs",
                         "gh_data/data/hairy3/odd_edges_odd_hairs"),
                     ("gh_data/data/hairy/odd_edges_even_hairs", "gh_data/data/hairy3/odd_edges_even_hairs")]:
    for file in get_g6_files(path2):
        fullfile1 = os.path.join(path1, file)
        fullfile2 = os.path.join(path2, file)
        if os.path.exists(fullfile1):
            with open(fullfile1) as f:
                txt1 = f.read()
            with open(fullfile2) as f:
                txt2 = f.read()
            print(
                f"Comparing {fullfile2}: {'OK' if txt1 == txt2 else 'ERROR*************************'}")
