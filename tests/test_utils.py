import xlranker as xlr

SMALL_PROTEIN = xlr.Protein(name="Small", abundance=1)
BIG_PROTEIN = xlr.Protein(name="Big", abundance=2)
MISSING_A = xlr.Protein(name="Missing A", abundance=None)
MISSING_B = xlr.Protein(name="Missing B", abundance=None)


def test_protein_order_with_one_null():
    # small protein should always come first
    assert xlr.set_protein_order(SMALL_PROTEIN, MISSING_A) == (SMALL_PROTEIN, MISSING_A)
    assert xlr.set_protein_order(MISSING_A, SMALL_PROTEIN) == (SMALL_PROTEIN, MISSING_A)


def test_protein_order_with_both_null():
    # Output should be same order as input
    assert xlr.set_protein_order(MISSING_B, MISSING_A) == (MISSING_B, MISSING_A)
    assert xlr.set_protein_order(MISSING_A, MISSING_B) == (MISSING_A, MISSING_B)


def test_protein_order_no_nulls():
    assert xlr.set_protein_order(SMALL_PROTEIN, BIG_PROTEIN) == (
        BIG_PROTEIN,
        SMALL_PROTEIN,
    )
    assert xlr.set_protein_order(BIG_PROTEIN, SMALL_PROTEIN) == (
        BIG_PROTEIN,
        SMALL_PROTEIN,
    )
    same_val_as_small = xlr.Protein("Same as Small", SMALL_PROTEIN.abundance)
    assert xlr.set_protein_order(same_val_as_small, SMALL_PROTEIN) == (
        same_val_as_small,
        SMALL_PROTEIN,
    )
    assert xlr.set_protein_order(SMALL_PROTEIN, same_val_as_small) == (
        SMALL_PROTEIN,
        same_val_as_small,
    )
