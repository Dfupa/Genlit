"""
GenLit modules package - Custom functions for biomedical literature mining.
"""

from .perplexity_search import perplexity_API_search
from .literature_fetch import search_and_fetch_pubmed, configure_entrez
from .clinvar_fetch import search_and_fetch_clinvar
from .flair_models import initialize_flair_models
from .hgnc_dedup import verify_genes_with_hgnc

__all__ = [
    "perplexity_API_search",
    "search_and_fetch_pubmed",
    "configure_entrez",
    "search_and_fetch_clinvar",
    "initialize_flair_models",
    "verify_genes_with_hgnc",
]

__version__ = "1.0.1"
