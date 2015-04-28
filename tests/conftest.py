import pytest


@pytest.fixture
def se_miso_ids():
    return ['chr1:100:200:+@chr1:300:400:+@chr1:500:600:+',
            'chr2:500:600:-@chr2:300:400:-@chr2:100:200:-']

@pytest.fixture
def mxe_miso_ids():
    return ['chr1:100:200:+@chr1:300:400:+@chr1:500:600:+@chr1:700:800:+',
            'chr2:700:800:-@chr2:500:600:-@chr2:300:400:-@chr2:100:200:-']

