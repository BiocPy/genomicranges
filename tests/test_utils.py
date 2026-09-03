import pytest
import numpy as np
from genomicranges.utils import sanitize_strand_vector, _sanitize_strand_search_ops, extract_groups_from_granges, _sanitize_vec
from genomicranges.GenomicRanges import GenomicRanges
from iranges import IRanges

def test_sanitize_strand_vector():
    with pytest.raises(ValueError):
        sanitize_strand_vector(None)

    with pytest.raises(ValueError):
        sanitize_strand_vector(np.array([[1, 2], [3, 4]]))

    with pytest.raises(ValueError):
        sanitize_strand_vector(np.array([2, 3]))

    with pytest.raises(ValueError):
        sanitize_strand_vector(["+", "a"])

    with pytest.raises(ValueError):
        sanitize_strand_vector([1, 2, 3])

    with pytest.raises(ValueError):
        sanitize_strand_vector([1.2, 3.4])

def test_sanitize_vec():
    masked = np.ma.masked_array([1, 2, 3], mask=[0, 1, 0])
    res = _sanitize_vec(masked)
    assert res == [1, None, 3]

def test_sanitize_strand_search_ops():
    # query_strand: +, -, *
    assert _sanitize_strand_search_ops("1", "1") == 1 # + + -> +
    assert _sanitize_strand_search_ops("1", "-1") is None # + - -> None
    assert _sanitize_strand_search_ops("1", "0") == 1 # + * -> +

    assert _sanitize_strand_search_ops("-1", "1") is None # - + -> None
    assert _sanitize_strand_search_ops("-1", "-1") == -1 # - - -> -
    assert _sanitize_strand_search_ops("-1", "0") == -1 # - * -> -

    assert _sanitize_strand_search_ops("0", "0") == 1 # * * -> +
    assert _sanitize_strand_search_ops("0", "-1") == -1 # * - -> -
    assert _sanitize_strand_search_ops("0", "1") is None

def test_extract_groups_from_granges():
    gr = GenomicRanges(seqnames=["chr1", "chr2", "chr1"], ranges=IRanges([1, 2, 3], [4, 5, 6]), strand=["+", "-", "+"])

    # ignore_strand=True
    groups = extract_groups_from_granges(gr, ignore_strand=True)
    assert len(groups) == 2
    assert groups[0][0] == "chr1"
    assert (groups[0][1] == np.array([0, 2])).all()

    # ignore_strand=False
    groups2 = extract_groups_from_granges(gr, ignore_strand=False)
    assert len(groups2) == 2
