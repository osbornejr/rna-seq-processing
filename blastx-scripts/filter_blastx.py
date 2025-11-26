import sys
import os
import csv

# -------------------------------
# Handle command line arguments
# -------------------------------
if len(sys.argv) < 2:
    print("Usage: python filter_blastx.py <blast_output> [output_csv]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) == 3 else os.path.splitext(input_file)[0] + "_summary.csv"

# -------------------------------
# Helper: choose best hit line
# -------------------------------

def select_hits(table_lines):
    prot_ids = []
    ## find first five Arabidopsis hits 
    #first, we need to find top of table
    line_counter = 0
    table_top = 0
    for line in table_lines:
        line_counter += 1
        if len(prot_ids) == 5:
            ##only need a max of 5 ids (?)
            return prot_ids
        if "Sequences producing significant alignments:" in line:
            table_top = line_counter
            continue
        if table_top>0 and line_counter > table_top:
            #check e-value
            split = line.split()
            if len(split) < 2:
                continue
            if len(split) >= 2 and float(split[-1]) > 0.05:
                #if evalue is too high, continueV
                continue
            #get species
            if "OS=" in line and "OX=" in line:
                species = line.split("OS=")[1].split("OX=")[0].strip()
            else:
                species = "Species undefined"
            if species == "Arabidopsis thaliana":
                prot_ids.append(split[0])
                continue
            continue
    if len(prot_ids)>0:
        #after checking all lines in table, if there is an AT hit, we are done.
        return prot_ids
    #otherwise get best of the rest? just AT for now
    return None

# -------------------------------
# Prepare CSV writer
# -------------------------------
out = open(output_file, "w", newline="")
writer = csv.writer(out)
writer.writerow([
    "QueryID","ProtID","Species"
])

# -------------------------------
# State variables
# -------------------------------
inblock = False
insubject = False
intable = False
table_lines = []
query_id = None
query_len = None
prot_ids = None
prot_table = []
evalue = None
species = None
prot_id = None
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
            insubject = False
            table_lines = []
            query_id = None
            query_len = None
            prot_ids = None
            prot_table = []
            evalue = None
            species = None
            prot_id = None

            parts = line.split()
            query_id = parts[1]
            query_len = parts[2].split("=")[1]
            continue

        if not inblock:
            continue

        # -------------------------------
        # First identify a case of no hits
        # -------------------------------
        if line.startswith("***** No hits found *****"):
        ## if at there were no hits found, we can stop here and write this as description 
            writer.writerow([
                query_id,
                "***** No hits found *****",
                ""
            ])
            # reset for next query
            inblock = False
            insubject = False
            table_lines = []
            query_id = None
            query_len = None
            prot_id = None
            prot_table = []
            evalue = None
            species = None
            prot_id = None
            continue
      
        # -------------------------------
        # Stop at each query >
        # -------------------------------
        if line.startswith(">"):
            insubject = True
            prot_id = line[1:].strip().split()[0]
            if "OS=" in line:
                species = line.split("OS=")[1].strip()
            continue
        
        if insubject:
            if "OS=" in line and "OX=" in line:
                species = line.split("OS=")[1].split("OX=")[0].strip()
                continue
            if "OS=" in line:
                species = line.split("OS=")[1].strip()
                continue
            if "OX=" in line:
                species_end = line.split("OX=")[0].strip()
                species_add = species + " " + species_end
                species = species_add.strip()
                continue
            if "Expect" in line:
                evalue = float(line.split("Expect")[1].split("=")[1].split(",")[0].strip())
            else:
                continue
            if evalue < 0.05:
                writer.writerow([
                    query_id,
                    prot_id,
                    species
                ])
                # reset for next query
                inblock = False
                insubject = False
                table_lines = []
                query_id = None
                query_len = None
                prot_id = None
                prot_table = []
                evalue = None
                species = None
                prot_id = None
                continue
            else:
                #reset for next subject
                insubject = False
                evalue = None
                species = None
                prot_id = None
                continue
        if "Effective search space used:" in line:
            writer.writerow([
                query_id,
                "***** No sig hits found *****",
                ""
            ])
            # reset for next query
            inblock = False
            insubject = False
            table_lines = []
            query_id = None
            query_len = None
            prot_id = None
            prot_table = []
            evalue = None
            species = None
            prot_id = None
            continue
                


            
       # # -------------------------------
       # # Whilst in table, get lines until first query >
       # # -------------------------------
       # if line.startswith(">"):
       #     # if line starts with >, table finished a few lines back
       #     prot_table = select_hits(table_lines)
       #     if prot_table:
       #         prot_ids = " ".join(prot_table)
       #         writer.writerow([
       #             query_id,
       #             query_len,
       #             prot_ids
       #         ])
       #         # reset for next query
       #         inblock = False
       #         table_lines = []
       #         query_id = None
       #         query_len = None
       #         prot_ids = None
       #         prot_table = []
       #         continue
       #     else:
       #         writer.writerow([
       #             query_id,
       #             query_len,
       #             "***** No sig hits found *****"
       #         ])
       #         # reset for next query
       #         inblock = False
       #         table_lines = []
       #         query_id = None
       #         query_len = None
       #         prot_ids = None
       #         prot_table = []
       #         continue

       # else:
       #     table_lines.append(line.strip())
       # continue


out.close()
print("CSV written to:", output_file)

