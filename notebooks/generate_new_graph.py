import os
import argparse
from graph import generate_random_graph


def main():
    """
    Generates new graphs
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="path to data directory",
                        type=str)
    parser.add_argument("-nn", help="number of nodes in new graph",
                        type=int)
    parser.add_argument("-od", help="number of OD pairs in new graph",
                        type=int)

    #TODO: manage edge cases when n_OD is too large for instance

    args = parser.parse_args()
    path = os.path.join('Data/', args.p)

    _, _ = generate_random_graph(
        args.nn, args.od, path)

    return


if __name__ == "__main__":
    main()
