#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
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
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import common
from nomad.datamodel.metainfo import public

m_package = Package(
    name='amber_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='amber.nomadmetainfo.json'))


class x_amber_mdin_method(MCategory):
    '''
    Parameters of mdin belonging to section method.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_amber_mdin_method'))


class x_amber_mdout_single_configuration_calculation(MCategory):
    '''
    Parameters of mdout belonging to section_single_configuration_calculation.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_amber_mdout_single_configuration_calculation'))


class x_amber_mdout_method(MCategory):
    '''
    Parameters of mdin belonging to section method.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_amber_mdout_method'))


class x_amber_mdout_run(MCategory):
    '''
    Parameters of mdin belonging to settings run.
    '''

    m_def = Category(
        categories=[public.settings_run],
        a_legacy=LegacyDefinition(name='x_amber_mdout_run'))


class x_amber_mdin_run(MCategory):
    '''
    Parameters of mdin belonging to settings run.
    '''

    m_def = Category(
        categories=[public.settings_run],
        a_legacy=LegacyDefinition(name='x_amber_mdin_run'))


class x_amber_section_input_output_files(MSection):
    '''
    Temperory variable to store input and output file keywords
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_amber_section_input_output_files'))


class x_amber_section_single_configuration_calculation(MSection):
    '''
    section for gathering values for MD steps
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_amber_section_single_configuration_calculation'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_amber_atom_positions_image_index = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        PBC image flag index.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_atom_positions_image_index'))

    x_amber_atom_positions_scaled = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        Position of the atoms in a scaled format [0, 1].
        ''',
        a_legacy=LegacyDefinition(name='x_amber_atom_positions_scaled'))

    x_amber_atom_positions_wrapped = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='meter',
        description='''
        Position of the atoms wrapped back to the periodic box.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_atom_positions_wrapped'))

    x_amber_dummy = Quantity(
        type=str,
        shape=[],
        description='''
        dummy
        ''',
        a_legacy=LegacyDefinition(name='x_amber_dummy'))

    x_amber_mdin_finline = Quantity(
        type=str,
        shape=[],
        description='''
        finline in mdin
        ''',
        a_legacy=LegacyDefinition(name='x_amber_mdin_finline'))

    x_amber_traj_timestep_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_timestep_store'))

    x_amber_traj_number_of_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_number_of_atoms_store'))

    x_amber_traj_box_bound_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_box_bound_store'))

    x_amber_traj_box_bounds_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_box_bounds_store'))

    x_amber_traj_variables_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_variables_store'))

    x_amber_traj_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_traj_atoms_store'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_amber_barostat_target_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        MD barostat target pressure.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_amber_barostat_target_pressure'))

    x_amber_barostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD barostat relaxation time.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_amber_barostat_tau'))

    x_amber_barostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD barostat type, valid values are defined in the barostat_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_amber_barostat_type'))

    x_amber_integrator_dt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD integration time step.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_integrator, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_integrator_dt'))

    x_amber_integrator_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD integrator type, valid values are defined in the integrator_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_integrator, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_integrator_type'))

    x_amber_periodicity_type = Quantity(
        type=str,
        shape=[],
        description='''
        Periodic boundary condition type in the sampling (non-PBC or PBC).
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_integrator, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_periodicity_type'))

    x_amber_langevin_gamma = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Langevin thermostat damping factor.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_langevin_gamma'))

    x_amber_number_of_steps_requested = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of requested MD integration time steps.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_integrator, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_number_of_steps_requested'))

    x_amber_thermostat_level = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat level (see wiki: single, multiple, regional).
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_thermostat_level'))

    x_amber_thermostat_target_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        MD thermostat target temperature.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_thermostat_target_temperature'))

    x_amber_thermostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD thermostat relaxation time.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_thermostat_tau'))

    x_amber_thermostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat type, valid values are defined in the thermostat_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_amber_thermostat_type'))


class section_atom_type(common.section_atom_type):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_atom_type'))

    x_amber_atom_type_element = Quantity(
        type=str,
        shape=[],
        description='''
        Element symbol of an atom type.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_atom_type_element'))

    x_amber_atom_type_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        van der Waals radius of an atom type.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_atom_type_radius'))


class section_interaction(common.section_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_interaction'))

    x_amber_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_interaction_atom_to_atom_type_ref'))

    x_amber_number_of_defined_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_amber_number_of_defined_pair_interactions'))

    x_amber_pair_interaction_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_amber_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_pair_interaction_atom_type_ref'))

    x_amber_pair_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['x_amber_number_of_defined_pair_interactions', 2],
        description='''
        Pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_pair_interaction_parameters'))


class section_molecule_interaction(common.section_molecule_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_molecule_interaction'))

    x_amber_molecule_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each molecule interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_molecule_interaction_atom_to_atom_type_ref'))

    x_amber_number_of_defined_molecule_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions within a molecule (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_amber_number_of_defined_molecule_pair_interactions'))

    x_amber_pair_molecule_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_defined_molecule_pair_interactions', 2],
        description='''
        Molecule pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_pair_molecule_interaction_parameters'))

    x_amber_pair_molecule_interaction_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_amber_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions within a molecule.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_pair_molecule_interaction_to_atom_type_ref'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_amber_program_version_date = Quantity(
        type=str,
        shape=[],
        description='''
        Program version date.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_version_date'))

    x_amber_xlo_xhi = Quantity(
        type=str,
        shape=[],
        description='''
        test
        ''',
        a_legacy=LegacyDefinition(name='x_amber_xlo_xhi'))

    x_amber_data_file_store = Quantity(
        type=str,
        shape=[],
        description='''
        Filename of data file
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_file_store'))

    x_amber_program_working_path = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_working_path'))

    x_amber_program_execution_host = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_execution_host'))

    x_amber_program_execution_path = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_execution_path'))

    x_amber_program_module = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_module'))

    x_amber_program_execution_date = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_execution_date'))

    x_amber_program_execution_time = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_program_execution_time'))

    x_amber_mdin_header = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_mdin_header'))

    x_amber_mdin_wt = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_mdin_wt'))

    x_amber_section_input_output_files = SubSection(
        sub_section=SectionProxy('x_amber_section_input_output_files'),
        repeats=False,
        a_legacy=LegacyDefinition(name='x_amber_section_input_output_files'))


class section_topology(common.section_topology):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_topology'))

    x_amber_input_units_store = Quantity(
        type=str,
        shape=[],
        description='''
        It determines the units of all quantities specified in the input script and data
        file, as well as quantities output to the screen, log file, and dump files.
        ''',
        a_legacy=LegacyDefinition(name='x_amber_input_units_store'))

    x_amber_data_bond_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_bond_types_store'))

    x_amber_data_bond_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_bond_count_store'))

    x_amber_data_angle_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_angle_count_store'))

    x_amber_data_atom_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_atom_types_store'))

    x_amber_data_dihedral_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_dihedral_count_store'))

    x_amber_data_angles_store = Quantity(
        type=str,
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_angles_store'))

    x_amber_data_angle_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_angle_list_store'))

    x_amber_data_bond_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_bond_list_store'))

    x_amber_data_dihedral_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_dihedral_list_store'))

    x_amber_data_dihedral_coeff_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_dihedral_coeff_list_store'))

    x_amber_masses_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_masses_store'))

    x_amber_data_topo_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_amber_data_topo_list_store'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_amber_section_single_configuration_calculation = SubSection(
        sub_section=SectionProxy('x_amber_section_single_configuration_calculation'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_amber_section_single_configuration_calculation'))


m_package.__init_metainfo__()
