import csv
from ete3 import Tree
import argparse


def get_args(argv=None):
    argParser = argparse.ArgumentParser()

    argParser.add_argument("-i", "--input-file", help="Input path to the tree file", required=True)
    argParser.add_argument("-o", "--output-path", help="Output path for new tree")
    argParser.add_argument("--bootstrap-table", help="Provide bootstrap value table", action="store_true")

    return argParser.parse_args(argv)


def read_tree(args):
    return Tree(args.input_file)


def assign_node_ids(node, counter, bootstrap_mapping):
    if node.is_leaf():
        return counter
    else:
        node.name = "InternalNode_" + str(counter).zfill(4)
        counter += 1
        for child in node.children:
            counter = assign_node_ids(child, counter, bootstrap_mapping)
        bootstrap_mapping[node.name] = node.support
        return counter


def write_output(args, tree):
    if args.output_path:
        output_file_tree_name = args.input_file.split("/")[-1].split(".")[0]
        output_file_tree = f"{args.args.output_path}{output_file_tree_name}.nwk"
    else:
        output_file_tree = f"{args.input_file.split('.')[0]}_AddIntNodes.nwk"

    tree.write(format=1, outfile=output_file_tree)
    print("Tree successfully exported to: ", output_file_tree)

    if args.bootstrap_table:
        output_file_table = f"{output_file_tree.split('.')[0]}_bootstraps.tsv"

        with open(output_file_table, mode="w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerow(["Internal_node_ID", "Bootstrap_value"])
            for node_id, bootstrap_value in bootstrap_mapping.items():
                writer.writerow([node_id, bootstrap_value])

        print("Table successfully exported to: ", output_file_table)


if __name__ == '__main__':
    args = get_args()

    tree = read_tree(args=args)

    bootstrap_mapping = {}
    assign_node_ids(tree, 1, bootstrap_mapping)

    write_output(args=args, tree=tree)
