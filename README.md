# GenLit 🧬📚
**GeneLit**: A Biomedical Literature Mining and Gene-Variant Curation Tool

<img align="right" src="/imgs/genlit_logo.png">

A Python tool for extracting, validating, and cross-checking disease-associated genes and variants from biomedical literature (PubMed/PMC) and ClinVar databases using Natural Language Processing with Named Entity Recognition (NER), with optional AI-powered cross-validation through the Perplexity API.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Advanced Configuration](#advanced-configuration)
- [Output Formats](#output-formats)
- [Dependencies](#dependencies)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Author](#author)

---

## Overview

GenLit automates the discovery and validation of disease-associated genes and genetic variants through:

1. **Literature Mining**: Searches PubMed and PMC for disease-relevant biomedical articles
2. **Entity Extraction**: Uses Flair's HunFlair2 biomedical NER model to identify genes and variants from abstracts and full text
3. **Variant Detection**: Extracts genetic variants in multiple HGVS formats (cDNA, protein, dbSNP)
4. **ClinVar Integration**: Cross-references findings with ClinVar pathogenic assertions
5. **AI Cross-Validation**: Optional Perplexity integration for literature-based fact-checking

Designed for bioinformaticians, researchers, and computational biologists conducting systematic gene discovery and validation workflows.

---

## Key Features

**Multi-Source Integration**
- PubMed/PMC full-text article retrieval
- ClinVar variant pathogenicity data
- Automated cross-referencing between sources

**Advanced NLP Processing using Flair**
- HunFlair2 biomedical NER model (BioNER)
- SciSpacy sentence and word tokenization
- Entity linking and normalization
- Gene and variant extraction with confidence scoring

**Flexible Variant Detection through Regular expression**
- Full HGVS format: `NM_003978.5(PSTPIP1):c.1034A>G (p.Tyr345Cys)`
- cDNA notation: `c.123A>G`, `c.123_124del`
- Protein changes: `p.Ala123Val`, `p.M1V`
- dbSNP IDs: `rs12345678`

**Perplexity AI Cross-Check** (Optional)
- Literature-based fact-checking of candidates
- Categorizes genes as: Confirmed, Uncertain, or Rejected
- Structured JSON output with curated verdicts
- Customizable temperature and token limits for reasoning control

**Flexible Output**
- CSV format for spreadsheet analysis
- TXT format with tabular layout
- Structured JSON for programmatic processing
- Console output with optional verbose logging

---

## Installation

### Prerequisites

- **Python**: 3.10 or higher
- **Conda**: For managing the environment ([install Miniconda](https://docs.conda.io/projects/miniconda/))
- **Internet Connection**: Required for PubMed/PMC API and Perplexity (optional)

### Step 1: Clone the Repository

```bash
git clone https://github.com/dfupa/GenLit.git
cd GenLit
```

### Step 2: Create Conda Environment

Using the provided `genlit.yml`:

```bash
conda env create -f genlit.yml
conda activate genelit
```

Alternatively, create a fresh environment with core dependencies:

```bash
conda create -n genelit python=3.10 -y
conda activate genelit
conda install -c bioconda -c conda-forge biopython flair scispacy pandas -y
pip install python-dotenv perplexityai
```

### Step 3: Install Biomedical Models

GenLit requires Flair's biomedical NER models. Download on first use (automatic) or manually:

```bash
python -c "
from flair.nn import Classifier
from flair.models import EntityMentionLinker

Classifier.load('hunflair2')
EntityMentionLinker.load('gene-linker')
EntityMentionLinker.load('disease-linker')
EntityMentionLinker.load('species-linker')
"
```

### Step 4: Configure API Keys

Create a `.env` file in the GenLit root directory:

```bash
echo "PERPLEXITY_API_KEY=your_api_key_here" > .env
```

Get your Perplexity API key from [Perplexity Labs](https://www.perplexity.ai/api/).

Alternatively, you can provide a path of .env with an API key. Check the --help section.

---

## Quick Start

### Basic Usage (No Perplexity)

```bash
python Genlit.py -q "Hidradenitis Suppurativa" -e "youremail@gmail.com"
```

Outputs results to console showing genes and variants from literature and ClinVar.

### With Perplexity Cross-Check with default parameters

```bash
python Genlit.py -q "Type 2 Diabetes" -e "youremail@gmail.com" --perplexity -v -o results.csv
```

Includes AI-powered validation and saves results to CSV.

### Search and fetch 100 papers

```bash
python Genlit.py -q "Crohn's Disease" -e "youremail@gmail.com" -r 100 -o results.txt
```

Outputs formatted table to `results.txt`.

---

## Usage

### Command-Line Interface

```
python Genlit.py -q "DISEASE_NAME" [OPTIONS]
```

### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--query` | `-q` | **Required.** Disease name (e.g., "Hidradenitis Suppurativa") |

### Optional Arguments

#### Output Options

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--output` | `-o` | STDOUT | Output file path (`.csv` or `.txt`). By default, STDOUT directly to the terminal |
| `--verbose` | `-v` | False | Verbose |

#### Literature Search Options

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--retmax` | `-r` | 10 | Max PubMed/PMC articles to fetch |
| `--email` | `-e` | None | Email for NCBI Entrez API (ie: youremail@gmail.com) |
| `--clinical-relevance` | `-c` | "Pathogenic" "Likely pathogenic" | ClinVar significance levels to include |

#### Perplexity Options

| Argument | Flag | Default | Description |
|----------|------|---------|-------------|
| `--perplexity` | `-p` | False | Enable AI cross-validation. Select to activate |
| `--temperature` | | 0.2 | Temperature parameter (0.0=deterministic, 1.0=creative) |
| `--max-tokens` | | 2048 | Max token length |
| `--model` | | sonar-pro | Model: "sonar" (fast) or "sonar-pro" (better reasoning) |
| `--env-path` | | .env | Path to custom `.env` file |

### Examples

#### Example 1: Basic Disease Search

```bash
python Genlit.py -q "Celiac Disease" -e "youremail@gmail.com"
```

**Output**: Lists genes and variants from recent literature and ClinVar.

#### Example 2: With Custom Clinical Relevance

```bash
python Genlit.py -q "Rheumatoid Arthritis" \
  --clinical-relevance "Pathogenic" "Likely pathogenic" "Uncertain significance" \
  -e "youremail@gmail.com" \
  -o ra_results.csv
```

**Output**: Includes variants with uncertain significance; saves to CSV.

#### Example 3: AI Cross-Check with Custom Temperature

```bash
python Genlit.py -q "Type 1 Diabetes" \
  -e "youremail@gmail.com" \
  --perplexity \
  --temperature 0.1 \
  --model sonar-pro \
  --verbose \
  -o t1d_validated.csv
```

**Output**: Low-temperature (deterministic) AI validation of genes; saves with Perplexity verdicts.

#### Example 4: Large-Scale Search with Different .env

```bash
python Genlit.py -q "Lupus" \
  -e "youremail@gmail.com" \
  --retmax 50 \
  --perplexity \
  --env-path /home/user/.env.production \
  -o lupus_comprehensive.csv
```

**Output**: 50 PubMed articles; uses production API keys.

#### Example 5: High-Determinism Curation

```bash
python Genlit.py -q "Marfan Syndrome" \
  -e "youremail@gmail.com" \
  --perplexity \
  --temperature 0.0 \
  --max-tokens 1500 \
  -v \
  -o Marfan_deterministic.txt
```

**Output**: Deterministic AI reasoning; verbose console output.

---

## Project Structure

```
GenLit/
├── Genlit.py                          # Main entry point
├── README.md                          # This file
├── genlit.yml                         # Conda environment specification
├── .env                               # API keys (create locally)
├── modules/                           # Custom module package
│   ├── __init__.py                    # Package initialization
│   ├── flair_models.py                # Flair NER model initialization
│   ├── literature_fetch.py            # PubMed/PMC retrieval
│   ├── clinvar_fetch.py               # ClinVar database queries
│   └── perplexity_search.py           # Perplexity API integration
├── imgs/                              # Images and logos
│   └── genlit_logo.png                # GenLit branding
└── .gitignore                         # Git ignore rules
```

### Core Modules

**`modules/flair_models.py`**
- Lazy-loads HunFlair2 NER, SciSpacy tokenizer, and entity linkers
- Configures CPU/GPU device support (CPU as default, edit if required)

**`modules/literature_fetch.py`**
- Searches NCBI PubMed/PMC via Entrez API
- Retrieves abstracts and full-text XML
- Handles Entrez API rate-limiting and error handling

**`modules/clinvar_fetch.py`**
- Queries ClinVar for disease-associated variants
- Filters by clinical significance
- Extracts gene-variant associations

**`modules/perplexity_search.py`**
- Interfaces with Perplexity Chat API
- Generates structured JSON verdicts
- Supports preset and control modes with customizable parameters

---

## Advanced Configuration

### Environment Variables

Create or modify `.env`:

```bash
# Required for Perplexity
PERPLEXITY_API_KEY=sk-...

# Optional NCBI Entrez settings. Feel free to proide the Entrez as a parameter " -e $ENTREZ_EMAIL "
ENTREZ_EMAIL=your.email@example.com
```

### Custom Temperature for Different Use Cases

- **0.0-0.2**: Best for systematic curation (deterministic)
- **0.3-0.5**: Balanced reasoning with some variation
- **0.6-1.0**: More creative; less reliable (not recommended for curation)

### Clinical Relevance Levels

ClinVar significance options:
- `Pathogenic`: Strong evidence of causation
- `Likely pathogenic`: Probable causation
- `Uncertain significance`: Unknown impact
- `Benign`: Not disease-causing
- `Likely benign`: Probably benign

And more: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/

---

## Output Formats

### CSV Output

Tab-separated values with columns:

| Column | Description |
|--------|-------------|
| Entity Type | Gene or Variant |
| Entity Name | Gene symbol or variant string |
| Associated Gene | Extracted gene (for variants) for genes is the same |
| Source | Abstract Only / ClinVar Only / Both / Full Text |
| Perplexity Status | Confirmed / Uncertain / Rejected / Not assessed |
| Notes | Additional notes |

**Example:**

```csv
Entity Type,Entity Name,Associated Gene,Source,Perplexity Status,Notes
Gene,PSTPIP1,PSTPIP1,Both,Confirmed,Found in literature (PubMed|PMC) + ClinVar
Variant,NM_003978.5(PSTPIP1):c.1034A>G,PSTPIP1,ClinVar Only,Confirmed,Found only in ClinVar but not in abstracts or full text
```

### TXT Output

Human-readable table format with header metadata:

```
Disease: Hidradenitis Suppurativa
Generated: 2026-01-28 14:30:45
================================================================================
Entity Type  Entity Name               Associated Gene  Source     Perplexity Status
Gene         NCSTN                    NCSTN            Both        Confirmed
Gene         PSTPIP1                  PSTPIP1          ClinVar...  Confirmed
Variant      NM_012275.3(IL36RN)...   IL36RN           Both        Confirmed
================================================================================
Total entities: 24
```

### Console Output (Verbose)

Real-time logging with timestamps:

```
[LOG] Initializing Flair models...
  ✓ SciSpacy sentence splitter loaded
  ✓ HunFlair2 NER tagger loaded
  ✓ Gene linker loaded
  ✓ Disease linker loaded
  ✓ Species linker loaded
[LOG] All models initialized successfully


================================================================================
[INFO] GeneLit Pipeline Started
================================================================================
[INFO] Disease Query: Hidradenitis Suppurativa
[INFO] PubMed retmax: 10
[INFO] Clinical relevance filter: Pathogenic, Likely pathogenic, Uncertain significance
[INFO] NCBI Email: youremail@gmail.com
[INFO] Output file: test.txt
[INFO] Perplexity cross-check: ENABLED
[INFO]   - Model: sonar-pro
[INFO]   - Temperature: 0.2
[INFO]   - Max tokens: 2048
[INFO] Env file: Default (.env)
================================================================================

Searching for: Hidradenitis Suppurativa



[LOG] Running Perplexity cross-check validation...

[LOG] Loaded .env from current directory


[INFO] Querying Perplexity for cross-validation...

[LOG] ✓ Successfully parsed structured JSON response
[LOG] Confirmed: 13 genes
[LOG] Uncertain: 3 genes
[LOG] Rejected: 23 genes
[LOG] Notes: Confirmed genes (NCSTN, PSENEN, PSTPIP1) and their listed variants have strong evidence from multiple sources, including MedlinePlus Genetics[1], comprehensive reviews[2], PubMed literature[3], and ClinVar listings for familial and sporadic HS cases involving gamma-secretase complex (NCSTN, PSENEN) and autoinflammatory pathways (PSTPIP1). KLF5 is uncertain due to emerging GWAS risk loci evidence from UNC studies[5][6], lacking ClinVar assertions or monogenic confirmation. KRT17 and GATA2 are uncertain: KRT17 variant in ClinVar but weak HS link (keratinization gene, no strong literature); GATA2 has indirect immunity links via consortium[8] but no direct HS association. All others rejected: lack credible peer-reviewed or ClinVar evidence for HS; represent markers (CD4, CD8, etc.), hormones (FSH, LH), unrelated gasdermins (GSDM family), or spurious text-mining hits without biological or genetic validation in HS literature[2][3][7]. Categorization prioritizes monogenic/familial evidence, GWAS loci, and ClinVar clinical assertions over unvalidated mentions.
```

---

## Dependencies

### Core Libraries

| Package | Version | Purpose |
|---------|---------|---------|
| `flair` | 0.15.1+ | Biomedical NER with HunFlair2 |
| `spacy` | 3.7.4 | NLP foundations library |
| `scispacy` | 0.6.2+ | Scientific text processing |
| `en_core_sci_sm` | 0.5.4+ | SciSpacy English models |
| `biopython` | 1.86+ | Bio.Entrez for PubMed API |
| `pandas` | 2.3.3+ | Big Data manipulation |
| `perplexityai` | 0.22.0+ | Perplexity API client |
| `python-dotenv` | 1.2.1+ | Environment variable management |

### Full Environment

See `genlit.yml` for the complete Conda environment specification with all transitive dependencies.

---

## Troubleshooting

### Issue: `ModuleNotFoundError: No module named 'modules'`

**Solution**: Ensure you're running `Genlit.py` from the GenLit root directory:

```bash
cd GenLit
python Genlit.py -q "Disease"
```

### Issue: Flair models not found

**Solution**: Download models manually:

```bash
python -c "from flair.nn import Classifier; Classifier.load('hunflair2')"
```

### Issue: `PERPLEXITY_API_KEY not found`

**Solution**: Create `.env` in GenLit root:

```bash
echo "PERPLEXITY_API_KEY=your_key_here" > .env
```

### Issue: PubMed API timeout

**Solution**: Reduce `--retmax` or increase NCBI rate limits:

```bash
python Genlit.py -q "Disease" -r 5  # Request fewer articles to begin with
```

### Issue: Out of memory with large results

**Solution**: Use CSV output (streaming) instead of console:

```bash
python Genlit.py -q "Disease" -o results.csv  # Lower memory footprint
```

### Issue: Perplexity API returns empty verdict

**Solution**: Check API key validity and request quota:

```bash
python Genlit.py -q "Disease" --perplexity --verbose  # See API errors and debug
```

---

## Performance Tips

- **For interactive use**: Keep `--retmax` ≤ 20 and skip Perplexity
- **For batch processing**: Set `--temperature 0.0` for maximum consistency (full deterministic)
- **For quick testing**: Use `--retmax 5` with `--model sonar` (cheaper)
- **For production**: Use `--temperature 0.1-0.2` with `--model sonar-pro`

**Typical runtimes:**
- Literature search: 90-180 seconds (10 articles)
- Variant extraction: 5-10 seconds
- Perplexity validation: 15-50 seconds
- **Total**: ~3-5 minutes per disease query in a normal 

---

## Citation

If you use GenLit in your research, please cite:

```bibtex
@software{genlit2026,
  author = {Fuentes Palacios, Diego},
  title = {GenLit: Biomedical Literature Mining and Gene-Variant Curation},
  year = {2026},
  url = {https://github.com/yourusername/GenLit}
}
```

---

## Author

**Diego Fuentes Palacios**  
Bioinformatician | Barcelona, Catalonia, Spain  
📧 [Email](mailto:diegofupa@gmail.com)  
🔗 [GitHub](https://github.com/dfupa)

---

## License

This project is licensed under the  GNU General Public License V3 - see LICENSE file for details.

---

## Acknowledgments

- **Flair** team for HunFlair2 biomedical NER model
- **NCBI** for PubMed/PMC API
- **ClinVar** for variant pathogenicity data
- **Perplexity AI** for reasoning API
- **SciSpacy** for scientific NLP tools

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## Support

For issues, questions, or feature requests:
- Open an [GitHub Issue](https://github.com/yourusername/GenLit/issues)
- Email: diegofupa@gmail.com

---

**Last Updated**: January 28, 2026  
**Version**: 1.0.0