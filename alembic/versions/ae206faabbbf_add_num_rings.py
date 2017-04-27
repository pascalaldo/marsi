"""add num rings

Revision ID: ae206faabbbf
Revises: ef39a4ae2c8c
Create Date: 2017-04-27 10:36:47.337419

"""
from IProgress import ProgressBar

from alembic import op
import sqlalchemy as sa

from marsi.io.db import Metabolite
from sqlalchemy.orm import sessionmaker

# revision identifiers, used by Alembic.
revision = 'ae206faabbbf'
down_revision = 'ef39a4ae2c8c'
branch_labels = None
depends_on = None

Session = sessionmaker()


def upgrade():
    op.add_column("metabolites", sa.Column("num_rings", sa.Integer, nullable=False))
    session = Session(op.get_bind())

    progress = ProgressBar()

    for metabolite in progress(session.query(Metabolite).yield_per(10000)):
        molecule = metabolite.molecule('openbabel', get3d=False)
        metabolite.num_rings = len(molecule.OBMol.GetLSSR())

    session.commit()


def downgrade():
    op.drop_column('metabolites', 'num_rings')
