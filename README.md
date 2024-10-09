# MG Extractor

MG Extractor is a tool for processing genomic data. The tool runs based on parameters in a configuration file and generates corresponding results.

## Directory Structure
```
.  
├── config  
├── example  
│   ├── assembly  
│   ├── cds  
│   ├── gene_lable2  
│   ├── hmm  
│   ├── protein  
│   └── species  
├── gutsmash_LSP  
│   ├── gene_lable2  
│   ├── results_GCF_000005845.2_ASM584v2_genomic  
│   ├── results_GCF_000010245.2_ASM1024v1_genomic  
│   ├── results_GCF_000019425.1_ASM1942v1_genomic  
│   ├── results_GCF_000025465.1_ASM2546v1_genomic  
│   ├── results_GCF_000026225.1_ASM2622v1_genomic  
│   ├── results_GCF_000164535.1_ASM16453v1_genomic  
│   ├── results_GCF_000191145.1_ASM19114v1_genomic  
│   ├── results_GCF_000195655.1_ASM19565v1_genomic  
│   ├── results_GCF_001997115.1_ASM199711v1_genomic  
│   ├── results_GCF_019331655.1_ASM1933165v1_genomic  
│   └── species  
├── gutsmash_TMA  
│   ├── gene_lable2  
│   ├── results_GCF_000005845.2_ASM584v2_genomic  
│   ├── results_GCF_000010245.2_ASM1024v1_genomic  
│   ├── results_GCF_000019425.1_ASM1942v1_genomic  
│   ├── results_GCF_000025465.1_ASM2546v1_genomic  
│   ├── results_GCF_000026225.1_ASM2622v1_genomic  
│   ├── results_GCF_000164535.1_ASM16453v1_genomic  
│   ├── results_GCF_000191145.1_ASM19114v1_genomic  
│   ├── results_GCF_000195655.1_ASM19565v1_genomic  
│   ├── results_GCF_001997115.1_ASM199711v1_genomic  
│   ├── results_GCF_019331655.1_ASM1933165v1_genomic  
│   └── species  
├── result_LSP  
│   ├── GCF_000005845.2  
│   ├── GCF_000010245.2  
│   ├── GCF_000019425.1  
│   ├── GCF_001997115.1  
│   └── hmm_output  
├── result_TMA  
│   ├── Core_gene_TMA  
│   ├── GCF_000005845.2  
│   ├── GCF_000010245.2  
│   ├── GCF_000019425.1  
│   ├── GCF_000191145.1  
│   ├── GCF_001997115.1  
│   ├── GCF_019331655.1  
│   └── hmm_output  
└── script  
```
## Running the Script

To run the MG Extractor script, use the following command:
```
python script/MG_extractor.py config/config.TMA.json  
or  
python script/MG_extractor.py config/config.LSP.json  
```
## Configuration File Explanation  

Below is an example of the config.json configuration file with explanations for each parameter:  
```
{  
  "genome_file": "/path/to/assembly",            ## Path to the genome file  
  "output_dir": "result",                        ## Output directory  
  "gene_lable2": "/path/to/gene_lable2",         ## Path to the gene label file  
  "gene_name": "TMA",                            ## Gene name  
  "cds_dir": "/path/to/cds",                     ## Directory for CDS files  
  "protein_dir": "/path/to/protein",             ## Directory for protein files  
  "species_path": "/path/to/species",            ## Path to the species file  
  "gutsmash_output": "gutsmash_output",          ## GutSMASH output directory  
  "gap": 5,                                      ## Gene cluster Gap parameter  
  "cpu": 16,                                     ## Number of CPUs to use  
  "e_value": 0.00001,                            ## E-value threshold  
  "domT": 0.5,                                   ## Domain threshold  
  "score": 200,                                  ## Score threshold  
  "hmm_files": [".hmm"],                         ## List of HMM files  
  "hmm_dir": "/path/to/hmm",                     ## Directory for HMM files    
  "gutsmash_path": "/path/to/run_gutsmash.py",   ## Path to the GutSMASH script  
  "seqkit_path": "/path/to/seqkit",              ## Path to the SeqKit tool  
  "hmmsearch_path": "/path/to/hmmsearch "        ## Path to the HMMsearch tool  
}  
```
## Parameter Explanations  
```
genome_file: Path to the genome file.  
output_dir: Output directory for storing result files.  
gene_lable2: Path to the gene label file.  
gene_name: Gene name used to identify specific genes.  
cds_dir: Directory for CDS (coding sequences) files.  
protein_dir: Directory for protein files.  
species_path: Path to the species file.  
gutsmash_output: GutSMASH output directory.  
gap: Gap parameter specifying the maximum interval between genes.  
cpu: Number of CPUs to use for parallel computation.  
e_value: E-value threshold for filtering alignment results.  
domT: Domain threshold for filtering HMM alignment results.  
score: Score threshold for filtering alignment results.  
hmm_files: List of HMM files to use.  
hmm_dir: Directory for storing HMM files.  
gutsmash_path: Path to the GutSMASH script.  
seqkit_path: Path to the SeqKit tool.  
hmmsearch_path: Path to the HMMsearch tool.  
```
## Example Configuration Files  

TMA Example Configuration File  
```
{  
  "genome_file": "example/assembly",    
  "output_dir": "result_TMA",    
  "gene_lable2": "example/gene_lable2",    
  "gene_name": "TMA",  
  "cds_dir": "example/cds",  
  "protein_dir": "example/protein",  
  "species_path": "example/species",  
  "gutsmash_output": "gutsmash_TMA",  
  "gap": 5,  
  "cpu": 16,  
  "e_value": 0.00001,  
  "domT": 0.5,  
  "score": 200,  
  "hmm_files": ["cntA.hmm","cntB.hmm","cutC.hmm","cutD.hmm"],  
  "hmm_dir": "/path/to/hmm",  
  "gutsmash_path": "/path/to/run_gutsmash.py",  
  "seqkit_path": "/path/to/seqkit",  
  "hmmsearch_path": "/usr/bin/hmmsearch"  
}  
```
LSP Example Configuration File  
```
{  
  "genome_file": "example/assembly",  
  "output_dir": "result_LSP",  
  "gene_lable2": "example/gene_lable2",  
  "gene_name": "LSP",  
  "cds_dir": "example/cds",  
  "protein_dir": "example/protein",  
  "species_path": "example/species",  
  "gutsmash_output": "gutsmash_LSP",  
  "gap": 5,  
  "cpu": 16,  
  "e_value": 0.00001,  
  "domT": 0.5,  
  "score": 200,  
  "hmm_files": ["cntA.hmm","cntB.hmm","cutC.hmm","cutD.hmm"],  
  "hmm_dir": "/path/to/hmm",  
  "gutsmash_path": "/path/to/run_gutsmash.py",  
  "seqkit_path": "/path/to/seqkit",  
  "hmmsearch_path": "/path/to/hmmsearch"  
}  
```
## License  
```
This project is licensed under the Apache License 2.0 - see the LICENSE file for details.  
```
