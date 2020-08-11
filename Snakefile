# ===== Pipeline for running cellranger ========================================


# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
import subprocess
import shutil
import glob
import os 
import re



# Parameters from config.yaml
RAW_DATA     = config["RAW_DATA"]
RESULTS      = config["RESULTS"]
RNA_SAMPLES  = config["RNA_SAMPLES"]
GROUPS       = config["GROUPS"]
GENOME       = config["GENOME"]
ADT_SAMPLES  = config["ADT_SAMPLES"]
ADT_REF      = config["ADT_REF"]
ANTIBODIES   = config["ANTIBODIES"]
VDJ_SAMPLES  = config["VDJ_SAMPLES"]
VDJ_REF      = config["VDJ_REF"]
MAX_JOBS     = config["MAX_JOBS"]
LSF_TEMPLATE = config["LSF_TEMPLATE"]



# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")



# Set sample/group names
FASTQ_REGEX = "_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz"

if RNA_SAMPLES:
    RNA_SAMPLES = [x.strip() for x in RNA_SAMPLES]
    SAMPLES     = RNA_SAMPLES

if ADT_SAMPLES:
    ADT_SAMPLES = [re.sub(" ", "", x) for x in ADT_SAMPLES]
    SAMPLES     = [re.sub(",", "_", x) for x in ADT_SAMPLES]
    ADT_REF     = _check_path(ADT_REF)

    if RNA_SAMPLES:
        SAMPLES = [x + "-" + y for x, y in zip(RNA_SAMPLES, SAMPLES)]

if VDJ_SAMPLES:
    VDJ_SAMPLES = [x.strip() for x in VDJ_SAMPLES]
    VDJ_REF     = _check_path(VDJ_REF)

if GROUPS:
    GROUPS = [re.sub(" ", "", x) for x in GROUPS]
    GROUPS = [re.sub(",", "-", x) for x in GROUPS]



# Check directory/file paths
RAW_DATA = [_check_path(x) for x in RAW_DATA]
RESULTS  = _check_path(RESULTS)
GENOME   = _check_path(GENOME)

if LSF_TEMPLATE:
    LSF_TEMPLATE = _check_path(LSF_TEMPLATE)
else:
    LSF_TEMPLATE = "lsf"

FASTQ_DIR = RESULTS + "/fastqs"

if not os.path.exists(FASTQ_DIR):
    os.makedirs(FASTQ_DIR)



# Final output files
rule all:
    input:
        expand(
            "{results}/logs/antibody_csv.out",
            results = RESULTS
        ),
        expand(
            "{results}/logs/{sample}_count.out", 
            results = RESULTS, sample = SAMPLES
        ),
        expand(
            "{results}/logs/{group}_aggr.out", 
            results = RESULTS, group = GROUPS
        ),
        expand(
            "{results}/logs/{vdj_sample}_vdj.out",
            results = RESULTS, vdj_sample = VDJ_SAMPLES
        ),
        expand(
            "{results}/count_metrics.csv",
            results = RESULTS
        )

include: "src/rules/cellranger.snake"



