"""add solubility

Revision ID: ef39a4ae2c8c
Revises: c1b678211e9d
Create Date: 2017-04-21 12:30:00.211355

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = 'ef39a4ae2c8c'
down_revision = 'c1b678211e9d'
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('metabolites', sa.Column('solubility', sa.Float))


def downgrade():
    op.drop_column('metabolites', 'solubility')
