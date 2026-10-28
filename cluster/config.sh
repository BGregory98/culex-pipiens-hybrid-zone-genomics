#!/bin/bash
# =============================================================================
# cluster/config.sh
# =============================================================================
# Paths, reference files, standalone binaries and conda environment names
# shared by the cluster scripts. A script that needs them sources this file
# before anything else:
#
#   source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1
#
# Jobs were submitted from the project root with this directory installed as
# <project root>/scripts/, so the default path resolves from there. To submit
# from elsewhere, export CLUSTER_CONFIG=/absolute/path/to/config.sh.
#
# The project root is site-specific and is not stored here. Set PROJECT_ROOT
# in the environment, or in config.local.sh beside this file (untracked):
#
#   PROJECT_ROOT=/path/to/project
#
# Every value below derives from it and expands to the path used in the
# published run. Paths that the scripts give relative to the project root
# (locations.txt, logs/, data/..., software/minimap2-2.30/...) are left as
# they were run and do not appear here.
# =============================================================================

_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_config_dir}/config.local.sh" ]]; then
  source "${_config_dir}/config.local.sh"
fi
unset _config_dir

: "${PROJECT_ROOT:?PROJECT_ROOT is not set (see cluster/config.sh)}"

# ---- Reference genome: Cx. pipiens idCulPipi1 -------------------------------
REF_DIR="${PROJECT_ROOT}/reference_genomes"
# Contig names in pipe form (ENA|OZ004311|OZ004311.1), as in the BAM headers.
REF_FASTA="${REF_DIR}/cxpip_genome/cxpip_genome.fasta"
# Same sequence with | replaced by _ in contig names (02_mpileup/sanitize_headers.sh).
REF_FASTA_SANITIZED="${REF_DIR}/cxpip_genome/cxpip_genome.sanitized.fa"

# ---- Conda ------------------------------------------------------------------
CONDA_SH="${PROJECT_ROOT}/miniconda3/etc/profile.d/conda.sh"
ENV_POOLSNP="poolsnp"   # envs/poolsnp.yml
ENV_CPL2A="cpl2a"       # envs/cpl2a.yml

# ---- Standalone software ----------------------------------------------------
SOFTWARE_DIR="${PROJECT_ROOT}/software"
BWA_MEM2="${SOFTWARE_DIR}/bwa-mem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2"   # BWA-MEM2 2.2.1
JAVA_HOME_DIR="${SOFTWARE_DIR}/java/jdk-24.0.2"
PICARD_JAR="${SOFTWARE_DIR}/picard/3.4.0/picard.jar"                    # Picard 3.4.0
BAYPASS="${SOFTWARE_DIR}/baypass_public-master/sources/g_baypass"       # BayPass 3.1
BAYPASS_UTILS="${SOFTWARE_DIR}/baypass_public-master/utils/baypass_utils.R"
DELLY="${SOFTWARE_DIR}/delly/delly-v1.7.3-linux-amd64"                  # DELLY 1.7.3

# samtools 1.21 build used by the alignment step (sort, index, flagstat). It
# sits in a sibling directory of the project root.
SAMTOOLS_1_21_DIR="${SAMTOOLS_1_21_DIR:-$(dirname "${PROJECT_ROOT}")/North_etal_2023_reanalysis/programs/samtools-1.21}"
