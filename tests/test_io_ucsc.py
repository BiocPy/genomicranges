import pytest
from unittest.mock import patch
import pandas as pd
from genomicranges.io.ucsc import access_gtf_ucsc, read_ucsc
from genomicranges.GenomicRanges import GenomicRanges

def test_access_gtf_ucsc():
    url = access_gtf_ucsc("hg19", type="refGene")
    assert url == "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/genes//hg19.refGene.gtf.gz"

    with pytest.raises(ValueError):
        access_gtf_ucsc("hg19", type="invalidType")


@patch("genomicranges.io.ucsc.parse_gtf")
def test_read_ucsc(mock_parse_gtf):
    # Mock the return of parse_gtf with a dummy dataframe
    mock_df = pd.DataFrame({
        "seqnames": ["chr1"],
        "starts": [100],
        "ends": [200],
        "strand": ["+"]
    })
    mock_parse_gtf.return_value = mock_df

    gr = read_ucsc("hg19", type="refGene")
    
    assert isinstance(gr, GenomicRanges)
    assert len(gr) == 1
    assert gr.get_seqnames()[0] == "chr1"
    
    # ensure it was called properly
    mock_parse_gtf.assert_called_once_with(
        "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/genes//hg19.refGene.gtf.gz", 
        compressed=True
    )
