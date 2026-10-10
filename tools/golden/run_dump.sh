#!/usr/bin/env bash
# Regenerates tests/data/golden_v1/ from the v1.1 baseline.
#
# Usage: tools/golden/run_dump.sh [baseline-ref]
#
# The baseline library is built from a temporary worktree of <baseline-ref>
# (default: v1.1.0-baseline) with the host C++ compiler, so the dump always
# reflects the original code and not the working tree.
set -euo pipefail

repo_root="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
baseline_ref="${1:-v1.1.0-baseline}"
out_dir="$repo_root/tests/data/golden_v1"
timeout_seconds="${RANDLIB_DUMP_TIMEOUT:-60}"
cxx="${CXX:-g++}"
jobs="$(nproc)"

work_dir="$(mktemp -d)"
trap 'git -C "$repo_root" worktree remove --force "$work_dir/src" >/dev/null 2>&1 || true; rm -rf "$work_dir"' EXIT

git -C "$repo_root" worktree add --detach "$work_dir/src" "$baseline_ref" >/dev/null
src="$work_dir/src"
mkdir -p "$work_dir/obj"

sources=$(grep -oE '[A-Za-z/]+\.cpp' "$src/RandLib.pro")
for file in $sources; do echo "$file"; done |
    xargs -P "$jobs" -I{} sh -c "$cxx -std=c++17 -O2 -I'$src' -c '$src/{}' -o '$work_dir/obj/'\$(echo {} | tr / _).o"
"$cxx" -std=c++17 -O2 -I"$src" "$repo_root/tools/golden/dump_v1.cpp" "$work_dir"/obj/*.o -o "$work_dir/dump_v1"

rm -rf "$out_dir"
mkdir -p "$out_dir/parts"
"$work_dir/dump_v1" --list > "$work_dir/cases.tsv"

run_case() {
    local index="$1"
    local id
    id="$(awk -F'\t' -v i="$index" '$1 == i {print $2}' "$work_dir/cases.tsv")"
    local part="$out_dir/parts/$(printf '%04d' "$index").csv"
    local log="$work_dir/case_$index.log"
    local skip=""
    local attempt
    for attempt in $(seq 1 20); do
        local status=0
        timeout "$timeout_seconds" "$work_dir/dump_v1" --case "$index" --skip "$skip" > "$part" 2> "$log" || status=$?
        [ "$status" -eq 0 ] && return
        if [ "$status" -ne 124 ]; then
            printf '"%s",__case__,,,crash:%s\n' "$id" "$status" >> "$part"
            return
        fi
        skip="$skip,$(grep '^CALL ' "$log" | tail -1 | cut -d' ' -f2-)"
    done
    printf '"%s",__case__,,,timeout\n' "$id" >> "$part"
}
export -f run_case
export out_dir work_dir timeout_seconds

cut -f1 "$work_dir/cases.tsv" | xargs -P "$jobs" -I{} bash -c 'run_case {}'

while IFS=$'\t' read -r index id; do
    class="${id%%(*}"
    target="$out_dir/$class.csv"
    [ -f "$target" ] || echo "case,function,arg,value,status" > "$target"
    cat "$out_dir/parts/$(printf '%04d' "$index").csv" >> "$target"
done < "$work_dir/cases.tsv"
rm -rf "$out_dir/parts"

echo "Wrote $(ls "$out_dir" | wc -l) files to $out_dir"
grep -h ',__case__,' "$out_dir"/*.csv || true
