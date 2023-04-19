import argparse
import json


def extract_insertions(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        insertions = node["muts"]
        if insertions:
            data[name] = insertions
    return data


def label_string(label_dict, include='all'):
    if include != 'all':
        label_dict = {k: v for k, v in label_dict.items() if include in k}
    # print(label_dict)
    label_string = "; ".join(
        [
            f'{k}: {", ".join(v)}'
            for k,v in label_dict.items()
        ]
    )
    return label_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add insertion branch labels",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input", type=str, metavar="JSON", required=False, help="input Auspice JSON")
    parser.add_argument("--insertions", type=str, required=False, help="insertions node data file")
    parser.add_argument(
        "--output", type=str, metavar="JSON", required=False, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        auspice_json = json.load(f)


    if args.insertions:
        with open(args.insertions, "r") as f:
            insertions = extract_insertions(json.load(f))
    else:
        insertions = {}

    def attach_labels(n):  # closure

        if n["name"] in insertions:
            if "branch_attrs" not in n:
                n["branch_attrs"] = {}
            if "labels" not in n["branch_attrs"]:
                n["branch_attrs"]["labels"] = {}
            if "aa" in n["branch_attrs"]["labels"]:
                n["branch_attrs"]["labels"]["aa"] += "; "
            else:
                n["branch_attrs"]["labels"]["aa"] = ""
            insertion_string = label_string(insertions[n["name"]], include='ins')
            all_string = label_string(insertions[n["name"]], include='frame')
            n["branch_attrs"]["labels"]["aa"] += all_string
            n["branch_attrs"]["mutations"].update(insertions[n["name"]])

            n["branch_attrs"]["labels"]["insertions"] = insertion_string
            if n["name"].startswith("NODE_"):
                n["branch_attrs"]["labels"]["insertions (internal)"] = insertion_string

        if "children" in n:
            for c in n["children"]:
                attach_labels(c)

    attach_labels(auspice_json["tree"])

    with open(args.output, "w") as f:
        json.dump(auspice_json, f, indent=2)
