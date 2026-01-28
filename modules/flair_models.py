"""Flair NER model initialization wrapper module."""

import warnings
import flair
from flair.models import EntityMentionLinker
from flair.tokenization import SciSpacyTokenizer
from flair.splitter import SciSpacySentenceSplitter
from flair.nn import Classifier


def initialize_flair_models(verbose: bool = False) -> tuple:
    """
    Lazy-load Flair NER models and linkers.
    The model in question is based on hunflair2, accessible from Flair.
    It is a biomedical NER tagger trained on the BioNLP13CG corpus
    and other datasets for better performance in biomedical text.

    Args:
        verbose: Print loading status messages if True

    Returns:
        Tuple of (splitter, tagger, gene_linker, disease_linker, species_linker)
    """
    if verbose:
        print("[LOG] Initializing Flair models...")

    warnings.filterwarnings("ignore", category=FutureWarning)
    flair.device = "cpu"  # Use CPU for compatibility

    try:
        # Initialize sentence splitter
        splitter = SciSpacySentenceSplitter()
        if verbose:
            print("  ✓ SciSpacy sentence splitter loaded")

        # Load biomedical NER tagger
        tagger = Classifier.load("hunflair2")
        if verbose:
            print("  ✓ HunFlair2 NER tagger loaded")

        # Load entity linkers for normalization
        gene_linker = EntityMentionLinker.load("gene-linker")
        if verbose:
            print("  ✓ Gene linker loaded")

        disease_linker = EntityMentionLinker.load("disease-linker")
        if verbose:
            print("  ✓ Disease linker loaded")

        species_linker = EntityMentionLinker.load("species-linker")
        if verbose:
            print("  ✓ Species linker loaded")

        if verbose:
            print("[LOG] All models initialized successfully\n")

        return splitter, tagger, gene_linker, disease_linker, species_linker

    except Exception as e:
        print(f"\n[ERROR ]✗ While loading Flair models: {str(e)}")
        print("  Make sure you've installed the required packages:")
        print("  pip install -U flair scispacy en_core_sci_sm spacy==3.7.4")
        print("  pip install -U scispacy")
        print("  pip install -U en_core_sci_sm")
        print("  pip install flair")
        raise SystemExit(1)
