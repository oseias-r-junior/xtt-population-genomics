#!/usr/bin/env bash
#
# install_padloc.sh
#
# Creates the PADLOC environment, downloads the database,
# validates the installation and performs a smoke test.
#
# Author:
#     Oseias Rodrigues Feitosa Junior
#

set -euo pipefail

################################################################################
# Configuration
################################################################################

ENV_NAME="xtt_padloc"

PROJECT="${PROJECT_DIR:-$PWD}"

ENV_FILE="$PROJECT/workflow/envs/xtt_padloc.yml"

DB_DIR="$PROJECT/xtt_project/defense_systems/padloc/database"

TEST_GENOME="$PROJECT/xtt_project/annotations/fna/AB22007.fna"

TEST_OUT="$PROJECT/xtt_project/defense_systems/padloc/test"

################################################################################
# Micromamba
################################################################################

eval "$(micromamba shell hook --shell bash)"

################################################################################
# Create environment
################################################################################

echo
echo "====================================================="
echo "Creating PADLOC environment"
echo "====================================================="

micromamba env create -f "$ENV_FILE"

################################################################################
# Activate
################################################################################

micromamba activate "$ENV_NAME"

################################################################################
# Basic validation
################################################################################

echo
echo "====================================================="
echo "Checking executables"
echo "====================================================="

for exe in padloc hmmscan prodigal blastn minced
do
    command -v "$exe" >/dev/null \
        || { echo "ERROR: $exe not found."; exit 1; }

    echo "OK  $exe"
done

################################################################################
# Download database
################################################################################

echo
echo "====================================================="
echo "Downloading PADLOC database"
echo "====================================================="

mkdir -p "$DB_DIR"

padloc \
    --db-update \


PADLOC_DATA=$(padloc --help | grep "default" | sed -E "s/.*default '(.*)'.*/\1/")

################################################################################
# Smoke test
################################################################################

echo
echo "====================================================="
echo "Running smoke test"
echo "====================================================="

mkdir -p "$TEST_OUT"

padloc \
    --fna "$TEST_GENOME" \
    --outdir "$TEST_OUT" \
    --data "$DB_DIR"

echo
echo "Smoke test completed."

################################################################################
# Done
################################################################################

echo
echo "====================================================="
echo "PADLOC installation completed successfully."
echo "====================================================="