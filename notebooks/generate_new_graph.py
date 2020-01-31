import os
import argparse
from graph import generate_random_graph


def main():
    """
    Generates new graphs
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path to data directory",
                        type=str)
    parser.add_argument("n_nodes", help="number of nodes in new graph",
                        type=int)
    parser.add_argument("n_edges", help="number of edges in new graph",
                        type=int)
    parser.add_argument("n_OD", help="number of OD pairs in new graph",
                        type=int)

    #TODO: manage edge cases when n_OD is too large for instance

    args = parser.parse_args()
    path = os.path.join('Data/', args.path)

    _, _ = generate_random_graph(
        args.n_nodes, args.n_edges, args.n_OD, path)

    return


if __name__ == "__main__":
    main()
