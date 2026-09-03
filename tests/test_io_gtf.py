import pandas as pd
from genomicranges.io.gtf import parse_gtf, read_gtf, _parse_all_attribute
import pytest
from unittest.mock import patch
import io
import gzip

def test_parse_all_attribute():
    row = {"group": 'gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene";'}
    res = _parse_all_attribute(row)
    assert res["gene_id"] == "ENSG00000223972.5"
    assert res["transcript_id"] == "ENST00000456328.2"
    assert res["gene_type"] == "transcribed_unprocessed_pseudogene"


def test_parse_gtf():
    mock_df = pd.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "source": ["havana", "havana"],
        "feature": ["gene", "transcript"],
        "starts": [11869, 11869],
        "ends": [14409, 14409],
        "score": [".", "."],
        "strand": ["+", "+"],
        "frame": [".", "."],
        "group": ['gene_id "ENSG0"; transcript_id "ENST0";', 'gene_id "ENSG0"; transcript_id "ENST0";']
    })

    with patch("pandas.read_csv") as mock_read_csv:
        mock_read_csv.return_value = mock_df

        df = parse_gtf("dummy.gtf", compressed=False)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "gene_id" in df.columns
        assert df["ends"].iloc[0] == 14408  # Because it subtracts 1 in parse_gtf


def test_read_gtf():
    mock_df = pd.DataFrame({
        "seqnames": ["chr1", "chr1"],
        "source": ["havana", "havana"],
        "feature": ["gene", "transcript"],
        "starts": [11869, 11869],
        "ends": [14409, 14409],
        "score": [".", "."],
        "strand": ["+", "+"],
        "frame": [".", "."],
        "group": ['gene_id "ENSG0"; transcript_id "ENST0";', 'gene_id "ENSG0"; transcript_id "ENST0";']
    })
    with patch("pandas.read_csv") as mock_read_csv:
        mock_read_csv.return_value = mock_df

        gr = read_gtf("dummy.gtf")

        assert len(gr) == 2
        assert gr.get_seqnames()[0] == "chr1"
        assert gr.get_mcols().shape[1] > 0
