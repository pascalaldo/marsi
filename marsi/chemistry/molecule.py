# Copyright 2017 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os

from marsi.chemistry import openbabel
from marsi.chemistry import rdkit


VALID_FP_FORMATS = openbabel.fps + rdkit.fps


class Molecule(object):
    """
    Object representing a molecule.

    """

    @classmethod
    def from_sdf(cls, path_or_str):
        """
        Builds a molecule from a file or string.

        Parameters
        ----------
        path_or_str : str
            The input file name or a string with the molecule description in SDF format.

        Returns
        -------
        marsi.chemistry.molecule.Molecule

        """
        from_file = os.path.isfile(path_or_str)
        ob_mol = openbabel.sdf_to_molecule(path_or_str, from_file)
        ob_mol._settitle("")
        rd_mol = rdkit.sdf_to_molecule(path_or_str, from_file)
        return cls(ob_mol, rd_mol)

    @classmethod
    def from_inchi(cls, inchi):
        """
        Builds a molecule from an InChI key.

        Parameters
        ----------
        inchi : str
            A valid InChI key.

        Returns
        -------
        marsi.chemistry.molecule.Molecule

        """
        ob_mol = openbabel.inchi_to_molecule(inchi)
        rd_mol = rdkit.inchi_to_molecule(inchi)
        return cls(ob_mol, rd_mol)

    @classmethod
    def from_mol(cls, path_or_str):
        """
        Builds a molecule from a file or string.

        Parameters
        ----------
        path_or_str : str
            The input file name or a string with the molecule description in MOL format.

        Returns
        -------
        marsi.chemistry.molecule.Molecule

        """
        from_file = os.path.isfile(path_or_str)
        ob_mol = openbabel.mol_to_molecule(path_or_str, from_file)
        rd_mol = rdkit.mol_to_molecule(path_or_str, from_file)
        return cls(ob_mol, rd_mol)

    def __init__(self, ob_mol, rd_mol):
        self._ob_mol = ob_mol
        self._rd_mol = rd_mol

    @property
    def inchi(self):
        return openbabel.mol_to_inchi(self._ob_mol)

    @property
    def num_atoms(self):
        return self._ob_mol.OBMol.NumAtoms()

    @property
    def num_bonds(self):
        return self._ob_mol.OBMol.NumBonds()

    @property
    def num_rings(self):
        return len(self._ob_mol.OBMol.GetLSSR())

    def fingerprint(self, fpformat='maccs', bits=None):
        if fpformat not in VALID_FP_FORMATS:
            raise ValueError("Fingerprint '%s' is not valid. Use of of %s" % (fpformat, ", ".join(VALID_FP_FORMATS)))

        if fpformat in openbabel.fps:
            fp = openbabel.fingerprint(self._ob_mol, fpformat)
            bits = openbabel.fp_bits.get(fpformat, max(fp.bits))
            return openbabel.fingerprint_to_bits(fp, bits=bits)
        else:
            fp = rdkit.fingerprint(self._rd_mol, fpformat)
            if bits is None:
                bits = fp.GetNumBits()

            return rdkit.fingerprint_to_bits(fp, bits=bits)

    def _repr_html_(self):
        self._ob_mol.removeh()
        representation = self._ob_mol._repr_html_() or openbabel.mol_to_svg(self._ob_mol)
        self._ob_mol.addh()
        return representation
