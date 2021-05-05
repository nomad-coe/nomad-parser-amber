#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from amberparser import AmberParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return AmberParser()


def test_basic(parser):
    archive = EntryArchive()
    parser.parse('tests/data/polyAT_vac_init_min.out', archive, None)

    assert archive.section_run[0].program_version == '14 SANDER 2014'
    sec_system = archive.section_run[0].section_system[0]
    assert np.shape(sec_system.atom_positions) == (638, 3)
    assert sec_system.atom_labels[7] == 'O'
    sec_scc = archive.section_run[0].section_single_configuration_calculation
    assert sec_scc[1].energy_total.magnitude == approx(-3.49178376e-16)
