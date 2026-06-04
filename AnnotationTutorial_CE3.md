# Step 1: Merging Functional Annotation Files to build clean_annotations.tsv
```
#make sure you are in the correct directory that contains the target files
cd <your_working_directory>

#create a clean header
grep "^#query" CE3_functional.emapper.annotations | sed 's/#//' > clean_annotations.tsv

# Append data skipping comments
grep -v "^#" CE3_functional.emapper.annotations >> clean_annotations.tsv

# Verify it looks right
wc -l clean_annotations.tsv
head -3 clean_annotations.tsv | column -t
```
Once the ```clean_annotations.tsv``` file is created, good to move onto the next step

# Step 2: Quick BLAST Verification Strategy
Swissprot was used to generate the ```functional.emapper.annotations``` file, so it is always good to verify this against another system

## Install required modules
```
module load lang/Anaconda3

module load bio/seqtk/1.3-GCC-11.3.0
```

If there is an issue loading the package you might have to fully install. A couple ways is to see first if the module is available using ```module avail``` and ```module spider <packagename>``` to see if there is a version available to download

## Identify which annotations to BLAST
```
# View the column structure first
head -5 clean_annotations.tsv | column -t

# High COG/eggNOG score (adjust score column index as needed)
awk -F'\t' '$8 > 100' clean_annotations.tsv > high_score_annotations.tsv

# Also pull anything with specific functional categories of interest
# e.g., secreted proteins, virulence, transport, metabolism
grep -iE "virulence|secreted|effector|transporter|\
hormone|reproductive|gonadotropin|estrogen|progesterone|\
androgen|steroid|endocrine|vitellogenin|neuropeptide|\
GnRH|APGWamide|sex.determin|gametogenesis|oogenesis|\
spermatogenesis|follicle|luteinizing|prolactin|\
cytochrome.P450|aromatase|hydroxysteroid|\
shell|biomineralization|calcification|\
immune|defensin|lysozyme|antimicrobial|\
adhesion|foot.protein|mucin" \
clean_annotations.tsv >> high_score_annotations.tsv

# Deduplicate
sort -u high_score_annotations.tsv > priority_annotations.tsv

wc -l priority_annotations.tsv

# Check the ID format in your protein file
grep ">" CE3_highconfidence.aa | head -5

# Check the ID format in your priority list
head -5 priority_ids.txt
```
The ```grep``` command listing all of the categories of interest can be altered as you see fit! Just make sure to document the exact code you used.

## Extract the protein sequences for the priority genes
```
# Extract gene/protein IDs from priority list
cut -f1 priority_annotations.tsv > priority_ids.txt

# Pull sequences from your predicted proteome
seqtk subseq CE3_highconfidence.aa \
    priority_ids.txt > priority_proteins.faa

# Check the ID format in your extracted proteins
grep ">" priority_proteins.faa | head -5

# Check the IDs in your priority list
head -5 priority_ids.txt
```

## Load and Run BLAST
```
module load bio/BLAST+/2.13.0-gompi-2022a

# Common locations for shared databases on clusters
ls /opt/databases/ 2>/dev/null
ls /shared/db/ 2>/dev/null
ls /db/ 2>/dev/null
ls /reference/ 2>/dev/null

# Specifically look for SwissProt or nr
find /opt /shared /db /reference -name "swissprot*" -o -name "uniprot*" -o -name "nr.pal" 2>/dev/null | head -10

# Sort by score (column 4) descending and take top 500
sort -t$'\t' -k4,4 -rn priority_annotations.tsv | head -500 > top500_annotations.tsv

# Extract those IDs
cut -f1 top500_annotations.tsv > top500_ids.txt

# Extract their sequences
seqtk subseq CE3_highconfidence.aa \
    top500_ids.txt > top500_proteins.faa

# Verify
grep -c ">" top500_proteins.faa

module load bio/BLAST+/2.13.0-gompi-2022a

blastp \
  -query top500_proteins.faa \
  -db swissprot \
  -remote \
  -evalue 1e-10 \
  -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -out top500_blast_results.txt
```

# Step 3: Annotation Verification
Make sure you are in an interactive node before running the code! For the first step of the code we are actually going to be switching coding languages and running a python script!

## Filter the SwissProt BLAST Results
First we need to create the python script! Other than some minor adjustments, below is the code that will be copy and pasted into a new file
```
#create the document
nano filter_swissprot.py
```
Copy and paste the following
```
import sys
import os
import csv

def process_repro_blast(input_file, output_file):
    """
    Processes BLAST outfmt 6 without pandas.
    Format: qseqid sseqid pident length evalue bitscore stitle
    """
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # Dictionary to store the best hit for each gene
    # Key: qseqid, Value: [bitscore, evalue, full_row_list]
    best_hits = {}

    print(f"Reading BLAST results from: {input_file}")
    
    try:
        with open(input_file, 'r') as f:
            # Using tab delimiter for outfmt 6
            reader = csv.reader(f, delimiter='\t')
            
            count = 0
            for row in reader:
                if not row: continue
                count += 1
                
                # Extract columns based on your -outfmt string
                # 0:qseqid, 1:sseqid, 2:pident, 3:length, 4:evalue, 5:bitscore, 6:stitle
                try:
                    qseqid = row[0]
                    evalue = float(row[4])
                    bitscore = float(row[5])
                except (ValueError, IndexError):
                    continue # Skip malformed lines

                # Logic: If we haven't seen this gene, or if this hit is better
                # Higher bitscore is better. If bitscores equal, lower evalue is better.
                if qseqid not in best_hits:
                    best_hits[qseqid] = [bitscore, evalue, row]
                else:
                    current_best_bit = best_hits[qseqid][0]
                    current_best_eva = best_hits[qseqid][1]
                    
                    if (bitscore > current_best_bit) or (bitscore == current_best_bit and evalue < current_best_eva):
                        best_hits[qseqid] = [bitscore, evalue, row]

        # Write the filtered results to a new TSV
        header = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle']
        
        with open(output_file, 'w', newline='') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            writer.writerow(header)
            for gene in sorted(best_hits.keys()):
                writer.writerow(best_hits[gene][2])

        print(f"--- Processing Complete ---")
        print(f"Total lines processed: {count}")
        print(f"Unique genes with hits: {len(best_hits)}")
        print(f"Filtered results saved to: {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Usage: python filter_swissprot.py reproductive_blast_results.txt reproductive_top_hits.tsv
    input_fn = "reproductive_blast_results.txt"
    output_fn = "reproductive_top_hits.tsv"
    
    if len(sys.argv) > 1:
        input_fn = sys.argv[1]
    if len(sys.argv) > 2:
        output_fn = sys.argv[2]
        
    process_repro_blast(input_fn, output_fn)
```
Now lets run the script!
```
#Load Anaconda so that you are able to run python
module load lang/Anaconda3

#Run the python code
python filter_swissprot.py reproductive_blast_results.txt reproductive_top_hits.tsv
```
The above will run the full python script that will take the Raw BLAST results and filter them based on the top hits per gene.

## Merge Annotations
Like above we will have to create another python script. Like above, unless there are adjustments to file names that need to be made, it can be copy and pasted in after doing ```nano merge_annotations.py``` to create the document.
```
import sys
import os
import csv

def merge_annotations(swissprot_file, eggnog_file, output_file):
    if not os.path.exists(swissprot_file):
        print(f"Error: {swissprot_file} not found.")
        return
    if not os.path.exists(eggnog_file):
        print(f"Error: {eggnog_file} not found.")
        return

    # 1. Read eggNOG annotations into a dictionary
    print(f"Reading eggNOG data from: {eggnog_file}")
    eggnog_dict = {}
    
    with open(eggnog_file, 'r') as f_egg:
        for line in f_egg:
            # Skip header information lines
            if line.startswith('##'):
                continue
            
            # Find the column indices from the header row
            if line.startswith('#query'):
                headers = line.strip().lstrip('#').split('\t')
                try:
                    q_idx = headers.index('query')
                    name_idx = headers.index('Preferred_name')
                    desc_idx = headers.index('Description')
                except ValueError:
                    # Fallbacks if exact headers aren't found
                    q_idx, name_idx, desc_idx = 0, 4, 7 
                continue
            
            # Process actual data lines
            parts = line.strip('\n').split('\t')
            if len(parts) > q_idx:
                query_id = parts[q_idx]
                egg_name = parts[name_idx] if len(parts) > name_idx else "N/A"
                egg_desc = parts[desc_idx] if len(parts) > desc_idx else "N/A"
                eggnog_dict[query_id] = {'name': egg_name, 'desc': egg_desc}

    # 2. Read SwissProt hits and merge
    print(f"Reading SwissProt data from: {swissprot_file}")
    merged_data = []
    
    with open(swissprot_file, 'r') as f_swiss:
        reader = csv.reader(f_swiss, delimiter='\t')
        header = next(reader) # Grab SwissProt headers
        
        # Add our new eggNOG columns to the header
        new_header = header + ['eggNOG_Preferred_Name', 'eggNOG_Description']
        merged_data.append(new_header)
        
        for row in reader:
            if not row: continue
            qseqid = row[0]
            
            # Check if this gene is also in our eggNOG dictionary
            if qseqid in eggnog_dict:
                egg_name = eggnog_dict[qseqid]['name']
                egg_desc = eggnog_dict[qseqid]['desc']
            else:
                egg_name = "No eggNOG match"
                egg_desc = "No eggNOG match"
                
            merged_row = row + [egg_name, egg_desc]
            merged_data.append(merged_row)

    # 3. Write the final merged file
    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerows(merged_data)

    print(f"--- Merge Complete ---")
    print(f"Total genes processed from SwissProt: {len(merged_data) - 1}")
    print(f"Final merged table saved to: {output_file}")

if __name__ == "__main__":
    # Default filenames, can be overridden by command line arguments
    swiss_in = "reproductive_top_hits.tsv"
    egg_in = "CE3_functional.emapper.annotations"  # CE3 EggNOG annotations file
    out_file = "CE3_master_annotation_table.tsv"
    
    if len(sys.argv) > 1: swiss_in = sys.argv[1]
    if len(sys.argv) > 2: egg_in = sys.argv[2]
    if len(sys.argv) > 3: out_file = sys.argv[3]
        
    merge_annotations(swiss_in, egg_in, out_file)
```
Once the document is created you can run the script! (but always double check the file names)
```
python merge_annotations.py reproductive_top_hits.tsv CE3_functional.emapper.annotations CE3_master_annotation_table.tsv
```
The file generated (```CE3_master_annotation_table.tsv```) can be downloaded if you please as it is fun to scroll through!

# Step 4: Gene Annotation Search
Continuing to use python, this code will search the newly merged table for keywords that we have inputted. Once again if the file is not already created, that will have to be done. For this use ```nano search_genes.py``` and the code just below this can be copy and pasted in:
```
import pandas as pd
import os

def search_gene_annotations(file_path, search_terms):
    """
    Searches through a TSV annotation file for specific biological keywords.
    """
    if not os.path.exists(file_path):
        print(f"Error: Could not find the file '{file_path}'.")
        print("Make sure this script is in the same folder as your .tsv file.")
        return

    print(f"--- Loading {file_path} ---")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    print(f"Successfully loaded {len(df)} gene records.")
    df_searchable = df.astype(str).apply(lambda x: ' '.join(x), axis=1)

    results_list = []

    print("\n--- Starting Search for Reproductive Hormones ---")
    for term in search_terms:
        # Case-insensitive search
        mask = df_searchable.str.contains(term, case=False, na=False)
        matches = df[mask].copy()
        
        if not matches.empty:
            matches['Search_Term'] = term
            results_list.append(matches)
            print(f"✅ Found {len(matches)} matches for '{term}'")
        else:
            print(f"❌ No matches for '{term}'")

    if not results_list:
        print("\nNo genes matched any of the target reproductive hormones.")
        return

    # Combine all results into one table
    final_results = pd.concat(results_list).drop_duplicates(subset=[df.columns[0]])

    if not final_results.empty:
        output_file = "hormone_targets_results.csv"
        final_results.to_csv(output_file, index=False)
        print(f"\n--- SUCCESS ---")
        print(f"Total unique hormone candidates found: {len(final_results)}")
        print(f"Results saved to: {output_file}")
        
        # Display the first few matches
        print("\nPreview of matches:")
        cols_to_show = [df.columns[0]] + [c for c in df.columns if 'name' in c.lower() or 'desc' in c.lower()][:2] + ['Search_Term']
        print(final_results[cols_to_show].head(15))

if __name__ == "__main__":
    # 1. Name of your file — CE3 merged annotation table
    MY_FILE = "CE3_master_annotation_table.tsv"
    
    # 2. Targeted Reproductive Hormones
    # Using regex-like strings to catch variations (e.g., 'oxytocin' and 'vasotocin')
    KEYWORDS = [
        "conopressin",
        "oxytocin",
        "vasopressin",
        "vasotocin",      # Fish/bird variation often annotated in databases
        "egg laying",     # Catches "Egg laying hormone"
        "ELH",
        "gonadotropin",   
        "GnRH",
        "Egg"
    ]
    
    search_gene_annotations(MY_FILE, KEYWORDS)
```
And the code can be run using ```python search_genes.py```

# Step 5: Conducting TBLASTN
We are done, for now, with python scripts! This next step is where the bait files come in. We are using those bait files to search through our genome to see if we can pick something up with similar order. Always make sure that you are operating in an interactive session.

## Load the modules
```
module load lang/Anaconda3
module load bio/BLAST+/2.13.0-gompi-2022a
```

## Run the Code!
Once those are loaded we can run the script! The following can be copy and pasted in directly into the terminal, as long as the files and file paths are correct. If changes need to be made I suggest pasting it into a notepad file, making changes there and then copy and pasting it directly to the terminal.
```
# --- PATHS ---
# Make sure these point to your actual genome and the bait file you just made!
GENOME="CE3_r4_final.fasta.masked"
BAIT="ELH_bait.fasta"
OUT_DIR="<your_output_directory>"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "Step 1: Making BLAST database from your genome..."
# We create a temporary DNA database from your genome
makeblastdb -in "$GENOME" -dbtype nucl -out CE3_genome_db

echo "Step 2: Running TBLASTN with Gastropod ELH bait..."
# tblastn searches your protein bait against the DNA genome
# We use a very relaxed E-value (10.0) because peptide hormones evolve fast and are tiny
tblastn \
    -query "$BAIT" \
    -db CE3_genome_db \
    -evalue 10.0 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
    -out CE3_ELH_genomic_hits.tsv \
    -num_threads 8

echo "Hunt Complete! Check CE3_ELH_genomic_hits.tsv"
```

## Extract the DNA sequence from the results from the TBLASTN
The second to last step is to finally get the DNA sequence out from our genome sequence. This is done by using the coordinates found in the previous bait filtering step.
```
module load bio/SAMtools/1.16.1-GCC-11.3.0

# First, ensure your genome is indexed (you only run this once)
samtools faidx CE3_r4_final.fasta.masked

# Next, extract the region — replace scaffold name and coordinates with your actual hit coordinates from the TBLASTN results
samtools faidx CE3_r4_final.fasta.masked <scaffold_name>:<start>-<end>
```

# Step 6: Convert to Protein and prep for SignalP
Havent figured this out yet so steps are pending!
