#!/bin/sh
# Docker entrypoint script.

# Wait until Postgres is ready
while ! pg_isready -q -h $DB_HOST -p 5432 -U $DB_USER
do
  echo "$(date) - waiting for database to start"
  sleep 2
done

echo "$(date) - Running marsi container"

export PGPASSWORD=$DB_PASSWORD
psql -h database -c 'create database "$DB_NAME";' -U $DB_USER
psql -h database -U $DB_USER $DB_NAME -f tests/fixtures/marsi-db-schema.sql
python3 bin/restore_db.py
psql -h database -d $DB_NAME -c 'SELECT COUNT(*) FROM metabolites;' -U $DB_USER

pytest tests/test_io_db.py

# pytest --disable-warnings tests/test_bigg_api.py || true
# pytest --disable-warnings tests/test_chemistry.py || true
# pytest --disable-warnings tests/test_cobra.py || true
# pytest --disable-warnings tests/test_data_retrieval.py || true
# pytest --disable-warnings tests/test_optmet_design.py || true
# pytest --disable-warnings tests/test_post_processing.py || true
# pytest --disable-warnings tests/test_targets.py || true
# pytest --disable-warnings tests/test_utils.py || true
