import Fasta_gz_filereader
import pandas as pd

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

def extract_gc_content(nucleotides):
    count = 0
    for nucleotide in nucleotides:
        if nucleotide == "G" or nucleotide == "C":
            count = count + 1

    return count

def extract_length(nucleotides):
    count = len(nucleotides)
    return count


def main():
    fasta_gz_file = 'C:/Users/Aden Le/Documents/Cenik Lab/Research Project/Mousegenome/Unprocessed Files/appris_mouse_v2_selected.fa.gz'
    genes_dict = Fasta_gz_filereader.parse_fasta_gz(fasta_gz_file)
    gene_analysis_dict = {}

    for gene, overview in genes_dict.items():
        utr5_length = len(overview.get_utr5_transcript())
        utr5_gc = extract_gc_content(overview.get_utr5_transcript())
        cds_length = len(overview.get_cds_transcript())
        cds_gc = extract_gc_content(overview.get_cds_transcript())
        utr3_length = len(overview.get_utr3_transcript())
        utr3_gc = extract_gc_content(overview.get_utr3_transcript())
        genedata = Gene_analysis_data(gene, utr5_length, utr5_gc, cds_length, cds_gc, utr3_length, utr3_gc)
        gene_analysis_dict[gene] = genedata

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

    # Create a DataFrame
    df = pd.DataFrame(data)

    csv_path = "mycsvfile.csv"
    df.to_csv(csv_path, index=False)


main()

def datavalidation():

    fasta_gz_file = 'C:/Users/Aden Le/Documents/Cenik Lab/Research Project/Mousegenome/Unprocessed Files/appris_mouse_v2_selected.fa.gz'
    genes_dict = Fasta_gz_filereader.parse_fasta_gz(fasta_gz_file)

    for values in genes_dict.values():
        myset = set()
        cds_transcript = values.get_cds_transcript()
        startcodon = cds_transcript[0:3]
        myset.add(startcodon)
    if myset == {'ATG'}:
        print("The data has been validated")
    else:
        print("There is an error within the data set!")

#datavalidation()


    


