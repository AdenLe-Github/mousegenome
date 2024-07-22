import Fasta_gz_filereader  # Import custom module for reading gzipped FASTA files
import pandas as pd  # Import pandas for data manipulation and analysis

# Define a class to store gene analysis data
class Gene_analysis_data:
    def __init__(self, gene, utr5_length, utr5_gc, cds_length, cds_gc, utr3_length, utr3_gc):
        self.gene = gene
        self.utr5_length = utr5_length
        self.utr5_gc = utr5_gc
        self.cds_length = cds_length
        self.cds_gc = cds_gc
        self.utr3_length = utr3_length
        self.utr3_gc = utr3_gc
        self.total_length = utr5_length + cds_length + utr3_length

    # Method to print an overview of the gene's details
    def gene_overview(self):
        print(f"Gene Name: {self.gene}")
        print(f"UTR5 Length: {self.utr5_length}")
        print(f"UTR5 GC Count: {self.utr5_gc}")
        print(f"CDS Length: {self.cds_length}")
        print(f"CDS GC Count: {self.cds_gc}")
        print(f"UTR3 Length: {self.utr3_length}")
        print(f"UTR3 GC Count: {self.utr3_gc}")
        print(f"Transcript Size: {self.total_length}") 
        return ""

# Function to count the number of G and C nucleotides in a sequence
def extract_gc_content(nucleotides):
    count = 0
    for nucleotide in nucleotides:
        if nucleotide == "G" or nucleotide == "C":
            count = count + 1
    return count

# Function to get the length of a nucleotide sequence
def extract_length(nucleotides):
    count = len(nucleotides)
    return count

# Main function to execute the gene analysis workflow
def main():
    # File path to the gzipped FASTA file
    fasta_gz_file = 'C:/Users/Aden Le/Documents/Cenik Lab/Research Project/Mousegenome/Raw_data/appris_mouse_v2_selected.fa.gz'
    
    # Parse the FASTA file and get a dictionary of Gene objects
    genes_dict = Fasta_gz_filereader.parse_fasta_gz(fasta_gz_file)
    gene_analysis_dict = {}

    # Iterate over each gene and calculate the lengths and GC content of UTR5, CDS, and UTR3 regions
    for gene, overview in genes_dict.items():
        utr5_length = len(overview.get_utr5_transcript())
        utr5_gc = extract_gc_content(overview.get_utr5_transcript())
        cds_length = len(overview.get_cds_transcript())
        cds_gc = extract_gc_content(overview.get_cds_transcript())
        utr3_length = len(overview.get_utr3_transcript())
        utr3_gc = extract_gc_content(overview.get_utr3_transcript())
        
        # Create a Gene_analysis_data object and store it in the dictionary
        genedata = Gene_analysis_data(gene, utr5_length, utr5_gc, cds_length, cds_gc, utr3_length, utr3_gc)
        gene_analysis_dict[gene] = genedata

    # Prepare the data for DataFrame creation
    data = []
    for gene in gene_analysis_dict.values():
        data.append({
            "Gene Name": gene.gene,
            "UTR5 Length": gene.utr5_length,
            "UTR5 GC Count": gene.utr5_gc,
            "CDS Length": gene.cds_length,
            "CDS GC Count": gene.cds_gc,
            "UTR3 Length": gene.utr3_length,
            "UTR3 GC Count": gene.utr3_gc,
            "Transcript Size": gene.total_length
        })

    # Create a DataFrame from the data and save it to a CSV file
    df = pd.DataFrame(data)
    csv_path = "mouse_stats.csv"
    df.to_csv(csv_path, index=False)

# Execute the main function
main()

# Function to validate the data by checking the start codon of each CDS transcript
def datavalidation():
    # File path to the gzipped FASTA file
    fasta_gz_file = 'C:/Users/Aden Le/Documents/Cenik Lab/Research Project/Mousegenome/Raw_data/appris_mouse_v2_selected.fa.gz'
    
    # Parse the FASTA file and get a dictionary of Gene objects
    genes_dict = Fasta_gz_filereader.parse_fasta_gz(fasta_gz_file)

    # Check if the start codon of each CDS transcript is 'ATG'
    myset = set()
    for values in genes_dict.values():
        cds_transcript = values.get_cds_transcript()
        startcodon = cds_transcript[0:3]
        myset.add(startcodon)
    
    # Validate the data based on the start codon
    if myset == {'ATG'}:
        print("The data has been validated")
    else:
        print("There is an error within the data set!")

# Uncomment the following line to run data validation
# datavalidation()