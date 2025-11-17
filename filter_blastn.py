import sys
import os
import csv

# -------------------------------
# Handle command line arguments
# -------------------------------
if len(sys.argv) < 2:
    print("Usage: python filter_blast.py <blast_output> [output_csv]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) == 3 else os.path.splitext(input_file)[0] + "_summary.csv"

# -------------------------------
# RNA type keywords
# -------------------------------
rna_types_list = ["mrna", "lncrna", "rrna", "trna", "snrna", "scrna", "ncrna"]

# -------------------------------
# Helper: choose best hit line
# -------------------------------
def select_hit(table_lines):
    for line in table_lines:
        if "Cicer arietinum" in line:
            ##gets the highest rated Cicer transcript if possible
            return line
    ## if no cicer, we need to find top of table
    line_counter = 0
    table_top = 0
    for line in table_lines:
        line_counter += 1
        if "Sequences producing significant alignments:" in line:
            table_top = line_counter
            if line_counter == table_top + 2:
                ##gets highest rated general transcript
                return line
    return None

# -------------------------------
# Prepare CSV writer
# -------------------------------
out = open(output_file, "w", newline="")
writer = csv.writer(out)
writer.writerow([
    "QueryID","QueryLen","SeqID","Description",
    "RNAType","SeqLen","Score","Evalue"
])

# -------------------------------
# State variables
# -------------------------------
inblock = False
table_lines = []
description_lines = []

query_id = None
query_len = None
seq_id = None
description = None
rna_type = ""
seq_len = None
score = ""
evalue = ""

prev_line_blank = False
collecting_description = False

# -------------------------------
# Parse BLAST file
# -------------------------------
with open(input_file) as f:
    for line in f:

        # -------------------------------
        # New Query block
        # -------------------------------
        if line.startswith("Query="):
            inblock = True
            table_lines = []
            description_lines = []
            description = None
            rna_type = ""
            seq_len = None
            score = ""
            evalue = ""
            seq_id = None
            prev_line_blank = False
            collecting_description = False

            parts = line.split()
            query_id = parts[1]
            query_len = parts[2].split("=")[1]
            continue

        if not inblock:
            continue

        # -------------------------------
        # Collect table lines until line beginning with >
        # -------------------------------
        if not collecting_description:
            if line.startswith("***** No hits found *****"):
            ## if at there were no hits found, we can stop here and write this as description 
                writer.writerow([
                    query_id,
                    query_len,
                    "",
                    "***** No hits found *****",
                    "",
                    "",
                    "",
                    ""
                ])
                # reset for next query
                inblock = False
                table_lines = []
                description_lines = []
                description = None
                seq_id = None
                collecting_description = False
                prev_line_blank = False
                continue
            if line.startswith(">"):
                # if line starts with >, table finished a few lines back
                selected = select_hit(table_lines)
                if selected:
                    seq_id = selected.split()[0]
                    table_parts = selected.strip().split()
                    if len(table_parts) >= 2:
                        score = table_parts[-2]
                        evalue = table_parts[-1]
                    collecting_description = True
                    parts = line[1:].strip().split(None, 1)
                    if parts[0] == seq_id:
                        description_lines = [parts[1]] if len(parts) > 1 else []
                        continue
                    else:
                        continue
                    #need to check if first > line matches selected
            else:
                table_lines.append(line.strip())
                prev_line_blank = False
            continue

        # -------------------------------
        # Collect description lines until Length=
        # -------------------------------
        if collecting_description and seq_id:
            if line.startswith(">"):
                parts = line[1:].strip().split(None, 1)
                if parts[0] == seq_id:
                    description_lines = [parts[1]] if len(parts) > 1 else []
                    continue
                else:
                    continue 
            ##only look for Length= if description has begun... and append extra lines if not Length=
            if len(description_lines)>0:
                stripped = line.strip()
                if stripped.startswith("Length="):
                    description = " ".join(description_lines)
                    lower_desc = description.lower()
                    rna_type = ""
                    for t in rna_types_list:
                        if t in lower_desc:
                            rna_type = t.lower()
                            break
                    seq_len = stripped.split("=")[1]
                    writer.writerow([
                        query_id,
                        query_len,
                        seq_id,
                        description if description else "",
                        rna_type,
                        seq_len if seq_len else "",
                        score,
                        evalue
                    ])
                    # reset for next query
                    inblock = False
                    table_lines = []
                    description_lines = []
                    description = None
                    seq_id = None
                    collecting_description = False
                    prev_line_blank = False
                    continue
                else:
                    description_lines.append(stripped)
            else:
                continue

out.close()
print("CSV written to:", output_file)

