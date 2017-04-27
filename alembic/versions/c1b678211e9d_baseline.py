"""baseline

Revision ID: c1b678211e9d
Revises: 
Create Date: 2017-04-12 10:21:55.050636

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'c1b678211e9d'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "synonyms",
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('synonym', sa.String(500), nullable=False),
        sa.Index('uq_synonyms', 'synonym', unique=True)
    )

    op.create_table(
        'references',
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('database', sa.String(100), nullable=False),
        sa.Column('accession', sa.String(100), nullable=False),
        sa.UniqueConstraint('database', 'accession', name='_database_accession_uc')
    )

    op.create_table(
        'metabolites',
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('inchi_key', sa.String(27), nullable=False),
        sa.Column('inchi', sa.String(5000), nullable=False),
        sa.Column('analog', sa.Boolean, default=False),
        sa.Column('formula', sa.String(500), nullable=False),
        sa.Column('num_atoms', sa.Integer, nullable=False),
        sa.Column('num_bonds', sa.Integer, nullable=False),
        sa.Column('sdf', sa.Text, nullable=True),
        sa.Index('uq_inchi_key', 'inchi_key', unique=True),
    )

    op.create_table(
        'metabolite_fingerprints',
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('metabolite_id', sa.Integer, sa.ForeignKey('metabolites.id')),
        sa.Column('fingerprint_type', sa.String(10), nullable=False),
        sa.Column('fingerprint', sa.String(2048), nullable=False),
        sa.UniqueConstraint('metabolite_id', 'fingerprint_type', name='_fp_type_uc')
    )

    op.create_table(
        'metabolite_references',
        sa.Column('metabolite_id', sa.Integer, sa.ForeignKey('metabolites.id')),
        sa.Column('reference_id', sa.Integer, sa.ForeignKey('references.id'))
    )

    op.create_table(
        'metabolite_synonyms',
        sa.Column('metabolite_id', sa.Integer, sa.ForeignKey('metabolites.id')),
        sa.Column('synonym_id', sa.Integer, sa.ForeignKey('synonyms.id'))
    )


def downgrade():
    op.drop_table("metabolite_fingerprints")
    op.drop_table("metabolite_references")
    op.drop_table("metabolite_synonyms")
    op.drop_table("synonyms")
    op.drop_table("references")
    op.drop_table("metabolites")

