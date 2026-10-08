#!/usr/bin/env bash
# Run the complete Otter/Stavøya gref4hsi workflow on a protected working copy.

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
PYTHON="$REPOSITORY/gref4hsi_venv/bin/python"
RUNNER="$SCRIPT_DIR/process_otter_stavoya.py"

DEFAULT_SURVEY_ROOT="/media/leo/NESP_1/NTNU/UHI/Malin/2024_24_Oct_Otter_Stavøya/24_10_24"
RAW_DIR="${OTTER_RAW_DIR:-$DEFAULT_SURVEY_ROOT/raw}"
MISSION_DIR="${OTTER_MISSION_DIR:-$DEFAULT_SURVEY_ROOT/gref4hsi_full}"
STATE_DIR="$MISSION_DIR/.pipeline_state"
LOG_DIR="$MISSION_DIR/logs"
H5_DIR="$MISSION_DIR/Input/H5"

MIN_FREE_GIB=70

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -x "$PYTHON" ]] || die "Project interpreter not found: $PYTHON"
[[ -f "$RUNNER" ]] || die "Pipeline runner not found: $RUNNER"
[[ -d "$RAW_DIR" ]] || die "Raw directory not found: $RAW_DIR"
[[ "$RAW_DIR" != "$MISSION_DIR" ]] || die "Raw and mission directories must differ"

mapfile -d '' SOURCE_FILES < <(find "$RAW_DIR" -maxdepth 1 -type f -name '*.h5' -print0 | sort -z)
SOURCE_COUNT="${#SOURCE_FILES[@]}"
(( SOURCE_COUNT > 0 )) || die "No HDF5 files found in $RAW_DIR"

AVAILABLE_KIB="$(df -Pk "$RAW_DIR" | awk 'NR == 2 {print $4}')"
REQUIRED_KIB=$((MIN_FREE_GIB * 1024 * 1024))
if (( AVAILABLE_KIB < REQUIRED_KIB )); then
    die "Less than ${MIN_FREE_GIB} GiB is available on the mission disk"
fi

mkdir -p "$H5_DIR" "$STATE_DIR" "$LOG_DIR"

echo "Repository: $REPOSITORY"
echo "Raw input:  $RAW_DIR"
echo "Mission:    $MISSION_DIR"
echo "HDF5 files: $SOURCE_COUNT"
echo "Free space: $((AVAILABLE_KIB / 1024 / 1024)) GiB"
echo
echo "Creating copy-on-write clones when supported, ordinary copies otherwise."

for source in "${SOURCE_FILES[@]}"; do
    destination="$H5_DIR/$(basename -- "$source")"
    if [[ ! -e "$destination" ]]; then
        cp --reflink=auto --preserve=timestamps -- "$source" "$destination"
    fi
done

DESTINATION_COUNT="$(find "$H5_DIR" -maxdepth 1 -type f -name '*.h5' | wc -l)"
[[ "$DESTINATION_COUNT" -eq "$SOURCE_COUNT" ]] || \
    die "Working-copy verification failed: expected $SOURCE_COUNT files, found $DESTINATION_COUNT"

run_stage() {
    local stage="$1"
    local marker="$STATE_DIR/$stage.done"
    local log="$LOG_DIR/$stage.log"

    if [[ -f "$marker" ]]; then
        echo "Skipping completed stage: $stage"
        return
    fi

    echo
    echo "===== Running stage: $stage ====="
    MPLCONFIGDIR="${TMPDIR:-/tmp}/gref4hsi-matplotlib" \
        "$PYTHON" "$RUNNER" "$stage" --mission "$MISSION_DIR" \
        2>&1 | tee "$log"
    touch "$marker"
    echo "===== Completed stage: $stage ====="
}

run_stage prepare
run_stage pose
run_stage model
run_stage georeference
run_stage orthorectify

echo
echo "Pipeline completed successfully."
echo "RGB composites: $MISSION_DIR/Output/GIS/RGBComposites"
echo "Footprints:     $MISSION_DIR/Output/GIS/FootPrints"
echo "Logs:           $LOG_DIR"
