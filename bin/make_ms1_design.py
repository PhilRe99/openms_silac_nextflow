import sys

def main():

    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: build_ms1_groups.py "
            "<original_design.tsv> <runs.json> <group_design.tsv>"
        )

    original, runs, group = sys.argv[1:]

if __name__ == "__main__":
    main()
