import pytest


@pytest.fixture
def se_miso_ids():
    return ['chr1:100:200:+@chr1:300:400:+@chr1:500:600:+',
            'chr1:1500:1600:-@chr1:1300:1400:-@chr1:1100:1200:-',
            # This last event should be excluded from analyses because the test
            # fasta file only has a "chr1"
            'chr2:500:600:-@chr2:300:400:-@chr2:100:200:-']


@pytest.fixture
def mxe_miso_ids():
    return ['chr1:100:200:+@chr1:300:400:+@chr1:500:600:+@chr1:700:800:+',
            'chr1:1700:1800:-@chr1:1500:1600:-@chr1:1300:1400:-@chr1:1100:1200:-',
            # This last event should be excluded from analyses because the test
            # fasta file only has a "chr1"
            'chr2:700:800:-@chr2:500:600:-@chr2:300:400:-@chr2:100:200:-']


@pytest.fixture(params=['SE', 'MXE'])
def miso_ids_splice_type(request, se_miso_ids, mxe_miso_ids):
    if request.param == 'SE':
        return request, se_miso_ids
    elif request.param == 'MXE':
        return request, mxe_miso_ids
