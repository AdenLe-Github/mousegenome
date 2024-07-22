import gzip

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

    def gene_overview(self):
        print(f"Gene Name: {self.gene_name}")
        print(f"Transcript Size: {self.transcript_size}")
        print(f"UTR5 Region: {self.utr5_region}")
        print(f"CDS Region: {self.cds_region}")
        print(f"UTR3 Region: {self.utr3_region}")
        print(f"Transcript: {self.transcript[:50]}...")  # Print first 50 characters of the transcript for brevity
        return ""
    
    def get_utr5_transcript(self):
        return self.utr5_transcript

    def get_cds_transcript(self):
        return self.cds_transcript

    def get_utr3_transcript(self):
        return self.utr3_transcript
    
    def get_transcript_size(self):
        return self.transcript_size_calculated

    def get_transcript(self):
        return self.transcript

def parse_fasta_gz(file_path):
    genes = {}
    with gzip.open(file_path, 'rt') as file:
        header = None
        transcript_lines = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header and transcript_lines:
                    transcript = ''.join(transcript_lines)
                    try:
                        parts = header[1:].split('|')
                        gene_name = parts[5]
                        #transcript_size = int(parts[6])
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

                        # Debug print to check transcript length
                        #print(f"Gene: {gene_name}, Transcript Length: {len(transcript)}")
                    except (IndexError, ValueError) as e:
                        print(f"Skipping invalid header or data: {header}")
                
                header = line
                transcript_lines = []
            else:
                transcript_lines.append(line)
        
        # Handle the last entry
        if header and transcript_lines:
            transcript = ''.join(transcript_lines)
            try:
                parts = header[1:].split('|')
                gene_name = parts[5]
                #transcript_size = int(parts[6])
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

                # Debug print to check transcript length
                #print(f"Gene: {gene_name}, Transcript Length: {len(transcript)}")
            except (IndexError, ValueError) as e:
                print(f"Skipping invalid header or data: {header}")

    return genes

def isolate_regions(utr5, cds, utr3, transcript):
    if utr5 != "NA":
        index = utr5.find("-")
        lowerbound = 0
        upperbound = int(utr5[index + 1:])
        utr5_transcript = transcript[lowerbound:upperbound]
        transcript = transcript[upperbound:]

    elif utr5 == "NA":
        utr5_transcript = ""

    index = cds.find("-")
    lowerbound = int(cds[:index])
    upperbound = int(cds[index + 1:]) 
    #print(f"lowerbound is {lowerbound}")
    #print(f"upperbound is {upperbound}")
    index = (upperbound - lowerbound) + 1
    cds_transcript = transcript[:index]
    transcript = transcript[index:]

    if utr3 != "NA":
        index = utr3.find("-")
        lowerbound = int(utr3[:index])
        upperbound = int(utr3[index + 1:]) 
        index = (upperbound - lowerbound) + 1
        utr3_transcript = transcript[:index]   
        
        
    elif utr3 == "NA":
        utr3_transcript = ""

    return utr5_transcript, cds_transcript, utr3_transcript

