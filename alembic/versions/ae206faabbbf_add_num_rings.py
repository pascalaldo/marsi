"""add num rings

Revision ID: ae206faabbbf
Revises: ef39a4ae2c8c
Create Date: 2017-04-27 10:36:47.337419

"""
from alembic import op, util
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
    op.add_column("metabolites", sa.Column("num_rings", sa.Integer))
    util.messaging.log.info("adding rings to existing table")

    session = Session(bind=op.get_bind())
    for i, metabolite in enumerate(session.query(Metabolite).yield_per(10000)):
        molecule = metabolite.molecule('openbabel', get3d=False)
        metabolite.num_rings = len(molecule.OBMol.GetSSSR())

        if i % 1000 == 0:
            util.messaging.log.info("updated %i entries" % i)
            session.flush()

    session.commit()

    op.alter_column("metabolites", "num_rings", nullable=False)


def downgrade():
    op.drop_column('metabolites', 'num_rings')
