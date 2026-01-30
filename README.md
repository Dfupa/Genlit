# GenLit 🧬📚
**GeneLit**: A Biomedical Literature Mining and Gene-Variant Curation Tool

<img align="right" src="/imgs/genlit_logo.png">

*Disclaimer: Logo generated through ChatGPT Image Generator 1*

A Python tool for extracting, validating, and cross-checking disease-associated genes and variants from biomedical literature (PubMed/PMC) and ClinVar databases using Natural Language Processing with Named Entity Recognition (NER), with optional AI-powered cross-validation through the Perplexity API.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Perplexity Advanced Configuration](#perplexity-advanced-configuration)
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
- Exponential backoff with configurable retry attempts
- Toggle between cost-conscious (single call) and performance-optimized (concurrent) modes
- **Concurrent mode | Chunking mode**: Automatically splits large candidate lists to optimize API usage
  - Parallel API calls for faster processing (3-5x speedup)
  - Exponential backoff with configurable retry attempts

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
git clone https://github.com/dfupa/Genlit.git
cd Genlit
```

### Step 2: Create Conda Environment

Using the provided `genlit.yml`. Be wary though, sometimes installing through conda YAML files might fail:

```bash
conda env create -f genlit.yml
conda activate genelit
```

Alternatively, create a fresh environment with core dependencies:

```bash
conda create -n genelit python=3.10 -y
conda activate genelit
conda install -c bioconda -c conda-forge biopython flair scispacy pandas -y
pip install python-dotenv perplexityai flair
#Download directly the SciSpacy English models 
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_sm-0.5.4.tar.gz

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
python Genlit.py -q "Hidradenitis Suppurativa" -e "your.email@example.com"
```

Outputs results to console showing genes and variants from literature and ClinVar.

### With Perplexity Cross-Check with default parameters

```bash
python Genlit.py -q "Type 2 Diabetes" -e "your.email@example.com" --perplexity -v -o results.csv
```

Includes AI-powered validation and saves results to CSV.

### Search and fetch 100 papers

```bash
python Genlit.py -q "Crohn's Disease" -e "your.email@example.com" -r 100 -o results.txt
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
| `--email` | `-e` | **Required.** Email for NCBI Entrez API (ie: your.email@example.com) |

### Optional Arguments

#### Output Options

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--output` | `-o` | STDOUT | Output file path (`.csv` or `.txt`). By default, STDOUT directly to the terminal |
| `--verbose` | `-v` | False | Verbose output (shows API calls and reasoning) |

#### Literature Search Options

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--retmax` | `-r` | 10 | Max PubMed/PMC articles to fetch |
| `--clinical-relevance` | `-c` | "Pathogenic" "Likely pathogenic" | ClinVar significance levels to include |
| `--score` | `-s` | 0.75 | NER tagging minimum score for labeling entities (higher = better accuracy) |

#### AI Cross-Validation (Basic Options)

| Argument | Flag | Default | Description |
|----------|------|---------|-------------|
| `--perplexity` | `-p` | False | Enable Perplexity API cross-check. Select to activate |
| `--temperature` | | 0.2 | Temperature parameter (0.0=deterministic, 1.0=creative). Recommended: 0.0-0.3 for curation |
| `--model` | | sonar-pro | Model: "sonar" (fast, cheaper) or "sonar-pro" (better reasoning, default) |
| `--env-path` | | .env | Path to custom `.env` file containing PERPLEXITY_API_KEY |
| `--max-tokens` | | 2048 | Max response tokens per API call/chunk |

#### Chunking & Concurrency (Advanced Perplexity Options)

| Argument | Flag | Default | Description |
|----------|------|---------|-------------|
| `--enable-chunking` | | True | Enable smart chunking for large candidate lists. Set to False for cost-conscious mode (single API call) |
| `--max-concurrent` | | 3 | Maximum concurrent API calls when chunking enabled. Increase for faster processing (e.g., 5), decrease if rate-limited |
| `--max-retries` | | 3 | Maximum retry attempts per chunk for transient errors (rate limits, timeouts, 503/502 errors) |
| `--token-limit` | | 4096 | Total token budget for chunking strategy. Increase for fewer chunks, decrease to save tokens |

### Examples

#### Example 1: Basic Disease Search

```bash
python Genlit.py -q "Celiac Disease" -e "your.email@example.com"
```

**Output**: Lists genes and variants from recent literature and ClinVar.

#### Example 2: With Custom Clinical Relevance

```bash
python Genlit.py -q "Rheumatoid Arthritis" \
  --clinical-relevance "Pathogenic" "Likely pathogenic" "Uncertain significance" \
  -e "your.email@example.com" \
  -o ra_results.csv
```

**Output**: Includes variants with uncertain significance; saves to CSV.

#### Example 3: AI Cross-Check with Custom Temperature

```bash
python Genlit.py -q "Type 1 Diabetes" \
  -e "your.email@example.com" \
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
  -e "your.email@example.com" \
  --retmax 50 \
  --perplexity \
  --env-path /home/user/.env.production \
  -o lupus_comprehensive.csv
```

**Output**: 50 PubMed articles; uses production API keys.

#### Example 5: High-Determinism Curation

```bash
python Genlit.py -q "Marfan Syndrome" \
  -e "your.email@example.com" \
  --perplexity \
  --temperature 0.0 \
  --max-tokens 1500 \
  -v \
  -o Marfan_deterministic.txt
```

**Output**: Deterministic AI reasoning; verbose console output.

---

## Perplexity Advanced Configuration

GenLit's Perplexity module now includes **smart chunking**, **concurrent API calls**, and **robust retry logic** for handling large candidate lists efficiently. Choose between cost-conscious and performance-optimized modes.

### Overview: Chunking vs. Single API Call

| Feature | Chunking Enabled (Default) | Chunking Disabled |
|---------|---------------------------|-------------------|
| **API Calls** | Multiple (parallelizable) | Single |
| **Speed** | 3-5x faster (concurrent) | Slower (~30-60s) |
| **Cost** | Higher (more API calls) | Lower (1 call) |
| **Token Budget** | Uses `--token-limit` effectively | May hit token limits on large candidates |
| **Best For** | Large candidates dataset (300+ items), time-sensitive | Small candidates datasets (<100 items), cost-sensitive |
| **Not Assessed Rate** | Lower (all items processed) | Higher (items dropped if limit exceeded) |

### Usage Examples: Chunking Control

#### Mode 1: Cost-Conscious (Single API Call)

**Best for**: Small-to-medium candidate lists (< 200 items), budget constraints

```bash
python Genlit.py -q "Hidroadenitis Suppurativa" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking False \
  -o pancreas_cheap.csv
```

**Characteristics**:
- Single API call (~$0.01)
- ~30-40 seconds processing
- Best for quick evaluations

#### Mode 2: Fast Processing (Concurrent Chunks)

**Best for**: Large datasets (300-1000+ items), time-sensitive workflows

```bash
python Genlit.py -q "Breast Cancer" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking True \
  --max-concurrent 5 \
  -o breast_fast.csv
```

**Characteristics**:
- ~8-10 parallel API calls (processed 5 at a time)
- ~60-90 seconds total (3-5x speedup vs. single call)
- Lower % of not assessed items
- Higher cost (~$0.08-0.15)

#### Mode 3: Robust Rate-Limit Handling

**Best for**: Production pipelines with API rate limits

```bash
python Genlit.py -q "Lung Cancer" \
  -e "your.email@example.com" \
  --perplexity \
  --max-retries 5 \
  --enable-chunking True \
  -o lung_robust.csv
```

**Characteristics**:
- Up to 5 retry attempts per chunk (vs. default 3)
- Better handling of 429 (rate limit), 503 (service unavailable), timeout errors
- Exponential backoff: 1s, 2s, 4s, 8s, 16s between retries
- More resilient to API hiccups

#### Mode 4: Custom Token Budget for Large Datasets

**Best for**: Very large candidate lists (1000+ items)

```bash
python Genlit.py -q "Type 2 Diabetes" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking True \
  --token-limit 6000 \
  -o t2d_large.csv
```

**Characteristics**:
- 6000 token budget (vs. default 4096)
- Fewer chunks needed (more items per call)
- ~5-8 API calls instead of 10+
- Trade-off: Longer response times per chunk


### Chunking Strategy Details

**How it works**:

1. **Token Estimation**: Each gene/variant is tokenized (~4 chars = 1 token)
2. **Fixed Overhead**: Prompt template uses a fixed set of tokens
3. **Dynamic Chunking**: Items packed into chunks respecting `--token-limit`
4. **Safety Margin**: ~10-20% buffer reserved for response tokens
5. **Parallelization**: Chunks processed concurrently (limited by `--max-concurrent`)
6. **Aggregation**: Results deduplicated and merged into final verdict

**Example token calculation**:
- Disease query: ~100 tokens
- Prompt template: ~250 tokens
- Gene list (100 genes × 20 tokens each): ~2000 tokens
- **Total**: ~2350 tokens (below default 4096 limit)
- **Chunks needed**: 1 chunk

If total exceeds 4096:
- Split into 2+ chunks
- Each chunk respects remaining token budget
- Chunks processed in parallel (up to `--max-concurrent`)

### Retry Logic (Exponential Backoff)

When transient errors occur:

```
Attempt 1: Immediate → Error
  ↓ Wait 1.0 + random(0-1) = ~0.5-1.5s
Attempt 2: Retry → Error
  ↓ Wait 2.0 + random(0-1) = ~2.0-3.0s
Attempt 3: Retry → Error
  ↓ Wait 4.0 + random(0-1) = ~4.0-5.0s
Attempt 4: Retry → Success ✓
```

**Retriable Errors** (automatic retry):
- `429`: Rate limit exceeded
- `503`: Service unavailable
- `502`: Bad gateway
- Timeout errors
- Connection resets

**Non-Retriable Errors** (fail immediately):
- `401`: Invalid API key
- `400`: Malformed request
- Token limit exceeded
- Invalid JSON schema response

### Production-Ready Configurations

#### For Research (Balanced)

```bash
python Genlit.py -q "Disease" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking True \
  --max-concurrent 3 \
  --max-retries 3 \
  --token-limit 4096 \
  --temperature 0.2 \
  -o disease_research.csv
```

**Cost**: ~$0.05-0.10  
**Time**: ~60-90s  
**Reliability**: High (3 retries, smart chunking)

#### For High-Throughput Screening

```bash
python Genlit.py -q "Disease" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking True \
  --max-concurrent 5 \
  --max-retries 5 \
  --token-limit 6000 \
  --temperature 0.2 \
  -o disease_hts.csv
```

**Cost**: ~$0.08-0.15  
**Time**: ~45-60s  
**Reliability**: Very high (5 retries, larger token budget)

#### For Cost-Sensitive Analysis

```bash
python Genlit.py -q "Disease" \
  -e "your.email@example.com" \
  --perplexity \
  --enable-chunking False \
  --max-retries 3 \
  -o disease_cheap.csv
```

**Cost**: ~$0.01-0.02  
**Time**: ~30-60s  
**Reliability**: Moderate (single call, may have "Not assessed" items)

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
│   └── perplexity_search.py           # Perplexity API integration (with chunking)
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
- Interfaces with Perplexity Chat API with async support
- **Smart chunking**: Splits large candidate lists based on token budgets
- **Concurrent processing**: Parallel API calls with configurable concurrency limits
- **Retry logic**: Exponential backoff for transient errors
- **Result aggregation**: Deduplicates and merges chunk verdicts
- Generates structured JSON verdicts
- Supports cost-conscious mode (chunking disabled) and performance mode (concurrent chunks)

---

## Advanced Configuration

### Environment Variables

Create or modify `.env`:

```bash
# Required for Perplexity
PERPLEXITY_API_KEY=sk-...

# Optional NCBI Entrez settings. Feel free to provide the Entrez as a parameter " -e $ENTREZ_EMAIL "
ENTREZ_EMAIL=your.email@example.com
```

### Custom Temperature for Different Use Cases

- **0.0**: Fully deterministic (best for systematic curation)
- **0.1-0.2**: Conservative with minimal variation (recommended as Default)
- **0.3-0.5**: Balanced reasoning with some exploration
- **0.6-1.0**: More creative; less reliable (not recommended for curation)

### Clinical Relevance Levels

ClinVar significance options:
- `Pathogenic`: Strong evidence of causation
- `Likely pathogenic`: Probable causation
- `Uncertain significance`: Unknown impact
- `Benign`: Not disease-causing
- `Likely benign`: Probably benign

See [ClinVar Documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) for more details.

---

## Output Formats

### CSV Output

Tab-separated values with columns:

| Column | Description |
|--------|-------------|
| Entity Type | Gene or Variant |
| Entity Name | Gene symbol or variant string |
| Associated Gene | Extracted gene (for variants); for genes is the same |
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
[INFO] NCBI Email: your.email@example.com
[INFO] Output file: test.txt
[INFO] Perplexity cross-check: ENABLED
[INFO]   - Model: sonar-pro
[INFO]   - Temperature: 0.2
[INFO]   - Max tokens: 2048
[INFO]   - Enable chunking: True
[INFO]   - Max concurrent API calls: 3
[INFO]   - Max retries per chunk: 3
[INFO]   - Token limit for chunking: 4096
[INFO]   - Env file: Default (.env)
================================================================================

Searching for: Hidradenitis Suppurativa... [LOG] Running Perplexity cross-check validation...

[LOG] Loaded .env from current directory

[MODE] Chunking ENABLED - Aggressive mode (concurrent multi-chunk)
[LOG] Chunking 150 candidates based on 4096 token limit...
[LOG] Created 2 chunk(s) for processing
[LOG] Launching 2 concurrent queries (max 3 simultaneous)...

[LOG] ✓ Chunk 1 processed: 45 confirmed, 12 rejected, 5 uncertain, 3 not assessed.
[LOG] ✓ Chunk 2 processed: 38 confirmed, 8 rejected, 4 uncertain, 0 not assessed.
[LOG] Aggregated results: 83 confirmed, 20 rejected, 9 uncertain, 3 not assessed

[INFO] Querying Perplexity for cross-validation...

[LOG] ✓ Successfully parsed structured JSON response
[LOG] Confirmed: 83 genes
[LOG] Uncertain: 9 genes
[LOG] Rejected: 20 genes
[LOG] Notes: Confirmed genes (NCSTN, PSENEN, PSTPIP1) and their listed variants have strong evidence from multiple sources, including MedlinePlus Genetics[1], comprehensive reviews[2], PubMed literature[3], and ClinVar listings for familial and sporadic HS cases...
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
| `pandas` | 2.3.3+ | Data manipulation |
| `perplexityai` | 0.22.0+ | Perplexity API client (async support required) |
| `python-dotenv` | 1.2.1+ | Environment variable management |

### Full Environment

See `genlit.yml` for the complete Conda environment specification with all transitive dependencies.

---

## Troubleshooting

### Issue: `ModuleNotFoundError: No module named 'modules'`

**Solution**: Ensure you're running `Genlit.py` from the GenLit root directory:

```bash
cd GenLit
python Genlit.py -q "Disease" -e "your.email@example.com"
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
python Genlit.py -q "Disease" -e "your.email@example.com" -r 5  # Request fewer articles to begin with
```

### Issue: Out of memory with large results

**Solution**: Use CSV output (streaming) instead of console:

```bash
python Genlit.py -q "Disease" -e "your.email@example.com" -o results.csv  # Lower memory footprint
```

### Issue: Perplexity API returns "Not assessed" for many items

**Solution**: The token limit was exceeded. Try one of:

1. **Increase token limit** (allows more items per chunk):
   ```bash
   python Genlit.py -q "Disease" -e "your.email@example.com" --perplexity --token-limit 6000
   ```

2. **Reduce max-tokens** (shorter responses, more chunking):
   ```bash
   python Genlit.py -q "Disease" -e "your.email@example.com" --perplexity --max-tokens 1024
   ```

3. **Disable chunking for single call** (may lose data if truly over limit):
   ```bash
   python Genlit.py -q "Disease" -e "your.email@example.com" --perplexity --enable-chunking False
   ```

### Issue: Rate limit errors (429) or service errors (503)

**Solution**: Increase max retries and reduce concurrent calls:

```bash
python Genlit.py -q "Disease" -e "your.email@example.com" \
  --perplexity \
  --max-retries 5 \
  --max-concurrent 2 \
  --verbose
```

### Issue: Perplexity API returns invalid JSON

**Solution**: Check API response format and ensure API key is valid:

```bash
python Genlit.py -q "Disease" -e "your.email@example.com" \
  --perplexity \
  --verbose  # See raw API responses
```

---

## Performance Tips

- **For interactive use**: Keep `--retmax` ≤ 20, skip Perplexity or use `--enable-chunking False`
- **For batch processing**: Set `--temperature 0.0` for maximum consistency (fully deterministic)
- **For quick testing**: Use `--retmax 5` with `--model sonar` and `--enable-chunking False` (cheapest)
- **For production**: Use `--temperature 0.1-0.2` with `--model sonar-pro` and `--enable-chunking True`
- **For very large datasets (1000+ items)**: Use `--token-limit 6000` or `8000` with `--max-concurrent 5`

**Typical runtimes** (with Perplexity):
- Literature search: 90-180 seconds (10-50 articles)
- Variant extraction: 5-10 seconds
- Perplexity validation (chunking enabled): 60-120 seconds (concurrent)
- Perplexity validation (chunking disabled): 30-60 seconds (single call)
- **Total**: ~3-5 minutes per disease query

---

## Citation

If you use GenLit in your research, please cite:

```bibtex
@software{genlit2026,
  author = {Fuentes Palacios, Diego},
  title = {GenLit: Biomedical Literature Mining and Gene-Variant Curation},
  year = {2026},
  url = {https://github.com/dfupa/Genlit}
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

This project is licensed under the GNU General Public License V3 - see LICENSE file for details.

---

## Acknowledgments

- **Flair** team for HunFlair2 biomedical NER model
- **NCBI** for PubMed/PMC Entrez API
- **ClinVar** for variant pathogenicity data
- **Perplexity AI** for reasoning AI API 
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
- Open a [GitHub Issue](https://github.com/dfupa/Genlit/issues)
