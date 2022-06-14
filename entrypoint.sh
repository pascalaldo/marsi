#!/bin/sh
# Docker entrypoint script.

# Wait until Postgres is ready
while ! pg_isready -q -h $DB_HOST -p 5432 -U $DB_USER
do
  echo "$(date) - waiting for database to start"
  sleep 2
done

echo "$(date) - Running marsi container"

while getopts c:t:d: flag
do
    case "${flag}" in
        c) do_db_creation=true;;
        t) do_tests=true;;
        d) do_db_init=true;;
        l) do_loop=true;;
    esac
done

if $do_db_creation
then
  export PGPASSWORD=$DB_PASSWORD
  psql -h database -c "create database \"$DB_NAME\";" -U $DB_USER
  psql -h database -U $DB_USER $DB_NAME -f tests/fixtures/marsi-db-schema.sql
  python3 bin/restore_db.py
  psql -h database -d $DB_NAME -c 'SELECT COUNT(*) FROM metabolites;' -U $DB_USER
fi

if $do_tests
then
  pytest tests || true
fi

if $do_db_init
then
  marsi db init
fi

while $do_loop
do
  echo "$(date) - Running marsi container loop"
  sleep 10
done

exit 0
