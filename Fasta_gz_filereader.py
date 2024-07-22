import gzip

# Define a class to represent a gene with various attributes and methods
class Gene:
    def __init__(self, gene_name, transcript_size, utr5_region, utr5_transcript, cds_region, cds_transcript, utr3_region, utr3_transcript, transcript):
        self.gene_name = gene_name
        self.transcript_size = transcript_size
        self.utr5_region = utr5_region
        self.utr5_transcript = utr5_transcript
        self.cds_region = cds_region
        self.cds_transcript = cds_transcript
        self.utr3_region = utr3_region
        self.utr3_transcript = utr3_transcript
        self.transcript = transcript
        self.transcript_size_calculated = len(transcript)

    # Method to print an overview of the gene's details
    def gene_overview(self):
        print(f"Gene Name: {self.gene_name}")
        print(f"Transcript Size: {self.transcript_size}")
        print(f"UTR5 Region: {self.utr5_region}")
        print(f"CDS Region: {self.cds_region}")
        print(f"UTR3 Region: {self.utr3_region}")
        print(f"Transcript: {self.transcript[:50]}...")  # Print first 50 characters of the transcript for brevity
        return ""
    
    # Getter method for the UTR5 transcript
    def get_utr5_transcript(self):
        return self.utr5_transcript

    # Getter method for the CDS transcript
    def get_cds_transcript(self):
        return self.cds_transcript

    # Getter method for the UTR3 transcript
    def get_utr3_transcript(self):
        return self.utr3_transcript
    
    # Getter method for the calculated transcript size
    def get_transcript_size(self):
        return self.transcript_size_calculated

    # Getter method for the full transcript
    def get_transcript(self):
        return self.transcript

# Function to parse a gzipped FASTA file and extract gene information
def parse_fasta_gz(file_path):
    genes = {}
    with gzip.open(file_path, 'rt') as file:
        header = None
        transcript_lines = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Process the previous entry when a new header is found
                if header and transcript_lines:
                    transcript = ''.join(transcript_lines)
                    try:
                        parts = header[1:].split('|')
                        gene_name = parts[5]
                        transcript_size = len(transcript)
                        
                        # Initialize regions as "NA" (Not Available)
                        utr5_region = "NA"
                        cds_region = "NA"
                        utr3_region = "NA"
                        
                        # Extract region information from the header
                        for part in parts[7:]:
                            if part.startswith("UTR5:"):
                                utr5_region = part.split(':')[1]
                            elif part.startswith("CDS:"):
                                cds_region = part.split(':')[1]
                            elif part.startswith("UTR3:"):
                                utr3_region = part.split(':')[1]
                        
                        # Isolate the region transcripts from the full transcript
                        utr5_transcript, cds_transcript, utr3_transcript = isolate_regions(utr5_region, cds_region, utr3_region, transcript)
                        
                        # Create a Gene object and add it to the dictionary
                        gene = Gene(gene_name, transcript_size, utr5_region, utr5_transcript, cds_region, cds_transcript, utr3_region, utr3_transcript, transcript)
                        genes[gene_name] = gene
                    except (IndexError, ValueError) as e:
                        print(f"Skipping invalid header or data: {header}")
                
                # Update header and reset transcript lines for the new entry
                header = line
                transcript_lines = []
            else:
                transcript_lines.append(line)
        
        # Handle the last entry in the file
        if header and transcript_lines:
            transcript = ''.join(transcript_lines)
            try:
                parts = header[1:].split('|')
                gene_name = parts[5]
                transcript_size = len(transcript)
                
                utr5_region = "NA"
                cds_region = "NA"
                utr3_region = "NA"
                
                for part in parts[7:]:
                    if part.startswith("UTR5:"):
                        utr5_region = part.split(':')[1]
                    elif part.startswith("CDS:"):
                        cds_region = part.split(':')[1]
                    elif part.startswith("UTR3:"):
                        utr3_region = part.split(':')[1]
                
                utr5_transcript, cds_transcript, utr3_transcript = isolate_regions(utr5_region, cds_region, utr3_region, transcript)
                gene = Gene(gene_name, transcript_size, utr5_region, utr5_transcript, cds_region, cds_transcript, utr3_region, utr3_transcript, transcript)
                genes[gene_name] = gene
            except (IndexError, ValueError) as e:
                print(f"Skipping invalid header or data: {header}")

    return genes

# Function to isolate UTR5, CDS, and UTR3 regions from the full transcript
def isolate_regions(utr5, cds, utr3, transcript):
    if utr5 != "NA":
        index = utr5.find("-")
        lowerbound = 0
        upperbound = int(utr5[index + 1:])
        utr5_transcript = transcript[lowerbound:upperbound]
        transcript = transcript[upperbound:]
    else:
        utr5_transcript = ""

    index = cds.find("-")
    lowerbound = int(cds[:index])
    upperbound = int(cds[index + 1:])
    index = (upperbound - lowerbound) + 1
    cds_transcript = transcript[:index]
    transcript = transcript[index:]

    if utr3 != "NA":
        index = utr3.find("-")
        lowerbound = int(utr3[:index])
        upperbound = int(utr3[index + 1:])
        index = (upperbound - lowerbound) + 1
        utr3_transcript = transcript[:index]
    else:
        utr3_transcript = ""

    return utr5_transcript, cds_transcript, utr3_transcript
