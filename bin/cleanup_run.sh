#!/usr/bin/env bash
# Remove disposable local execution state; preserve published results and inputs.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bin/cleanup_run.sh [--delete] [--sudo]

With no flags, preview the cleanup without deleting anything.
  --delete  Delete local Nextflow work/cache, root Nextflow logs, and Python caches.
  --sudo    Use sudo for deletion of Docker-owned files; requires --delete.
  --help    Show this help.

Results (including their logs), RAW/FASTA inputs, metadata, source and historical
context are preserved. Deleting work/ and .nextflow/ prevents resuming old runs.
Stop the pipeline before deleting. This script requires Bash and Linux /proc
(including WSL) and resolves the repository relative to its own location.
EOF
}

die() {
    printf 'Cleanup refused: %s\n' "$*" >&2
    exit 1
}

delete=false
use_sudo=false
for arg in "$@"; do
    case "$arg" in
        --delete) delete=true ;;
        --sudo) use_sudo=true ;;
        --help|-h) usage; exit 0 ;;
        *) die "Unknown argument: $arg. Use --help." ;;
    esac
done
if "$use_sudo" && ! "$delete"; then
    die '--sudo requires --delete.'
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(cd -- "$script_dir/.." && pwd -P)
[[ "$repo_dir" != / && -f "$repo_dir/main.nf" && -f "$repo_dir/nextflow.config" ]] \
    || die 'The script must live in bin/ of this Nextflow repository.'

shopt -s nullglob
candidates=(
    "$repo_dir/work"
    "$repo_dir/.nextflow"
    "$repo_dir"/.nextflow.log*
    "$repo_dir/__pycache__"
    "$repo_dir/bin/__pycache__"
    "$repo_dir/tests/__pycache__"
)
targets=()
for path in "${candidates[@]}"; do
    if [[ -e "$path" || -L "$path" ]]; then
        targets+=("$path")
    fi
done

printf 'Repository: %s\n' "$repo_dir"
if [[ ${#targets[@]} == 0 ]]; then
    printf 'Nothing to clean.\n'
    exit 0
fi
printf 'Cleanup targets:\n'
for path in "${targets[@]}"; do
    # du does not follow symlinks. Its output is informational only.
    du -sh -- "$path" 2>/dev/null || printf '  %s (size unavailable)\n' "$path"
done
printf 'Published results and input data will be kept. Old runs will not be resumable.\n'
if ! "$delete"; then
    printf 'Preview only. Run with --delete to clean, or --delete --sudo for Docker-owned files.\n'
    exit 0
fi

# Inspect processes directly instead of matching a ps/grep command against itself.
# Block this repository's Nextflow JVM and any process still using its work tree,
# including task wrappers left behind after an interrupted Nextflow launcher.
[[ -r /proc/self/cmdline ]] || die 'Process checks require Linux /proc (or WSL).'
for proc in /proc/[0-9]*; do
    pid=${proc##*/}
    [[ "$pid" != "$$" && -r "$proc/cmdline" ]] || continue
    process_cwd=$(readlink -- "$proc/cwd" 2>/dev/null) || continue
    process_cwd=${process_cwd% (deleted)}
    if [[ "$process_cwd" == "$repo_dir/work" || "$process_cwd" == "$repo_dir/work/"* ]]; then
        die "PID $pid is still using work/. Stop its task or leave that directory first."
    fi
    process_args=()
    mapfile -d '' -t process_args < "$proc/cmdline" 2>/dev/null || continue
    command_line="${process_args[*]}"
    if [[ "$command_line" == *nextflow.cli.Launcher* || "$command_line" == *nextflow-*-one.jar* ]]; then
        if [[ "$process_cwd" == "$repo_dir" || "$command_line" == *"$repo_dir"* ]]; then
            die "Nextflow PID $pid is still active for this repository. Stop it first."
        fi
    fi
done

if "$use_sudo"; then
    command -v sudo >/dev/null || die 'sudo is not available.'
    sudo -v || die 'sudo authorization failed; nothing was deleted.'
else
    # Check directory permissions before starting, avoiding a partial cleanup
    # when a container has created unwritable prepared/ or qpx/ directories.
    for path in "${targets[@]}"; do
        [[ -d "$path" && ! -L "$path" ]] || continue
        if ! blocked=$(find "$path" -type d ! -writable -print -quit 2>/dev/null); then
            die 'Cannot inspect all directories. Rerun with --delete --sudo.'
        fi
        [[ -z "$blocked" ]] || die "Unwritable directory: $blocked. Rerun with --delete --sudo."
    done
fi

for path in "${targets[@]}"; do
    if "$use_sudo"; then
        sudo -- rm -rf -- "$path"
    else
        rm -rf -- "$path"
    fi
done
printf 'Cleanup complete. Results and inputs were preserved.\n'
