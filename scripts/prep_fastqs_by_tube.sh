#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: prep_fastqs_by_tube.sh -i sample_list.txt -o OUTDIR -m {symlink|concat}
                              [--report] [--strict] [-v|--verbose] [--dry-run]
                              [--skip-complete] [--force]

Groups FASTQs by dna_tube_id (supports '_d<digits>_' and '_pos<digits>_') and either:
  - symlink : Symlinks all lane FASTQs (R1/R2) into per-tube folders
  - concat  : Concatenates all lanes per read into tube_R1.fastq.gz / tube_R2.fastq.gz

Options:
  -i  Input file with FASTQ paths (one per line), e.g. sample_list_raw.txt (R1-only is OK)
  -o  Output base directory
  -m  Mode: 'symlink' or 'concat'
  --report       Write per-tube summary TSV/CSV under OUTDIR/report/
  --strict       Exit non-zero if any tube has mismatched lane counts or missing R1/R2
  -v|--verbose   Verbose logging (show parsing, grouping, per-tube actions)
  --dry-run      Show what would be done, without creating or modifying files
  --skip-complete Skip tubes that appear already completed (manifest + timestamps match)
  --force        Rebuild tubes even if they appear completed (overrides --skip-complete)

Notes:
  * If list contains only R1 paths, the script auto-adds matching R2 (if found on disk)
  * gzip concatenation is safe: 'cat lane1.fastq.gz lane2.fastq.gz > merged.fastq.gz'
  * This version hard-prevents duplicate inputs from being merged or manifested.
EOF
}

INPUT="" ; OUTDIR="" ; MODE=""
REPORT=0 ; STRICT=0 ; VERBOSE=0 ; DRYRUN=0 ; SKIP_COMPLETE=0 ; FORCE=0

log()   { echo -e "[$(date +'%H:%M:%S')] $*"; }
vlog()  { [[ "${VERBOSE}" -eq 1 ]] && log "$*"; }
doit()  { if [[ "${DRYRUN}" -eq 1 ]]; then echo "DRY-RUN: $*"; else eval "$@"; fi; }

# ---------- Parse args ----------
while (( "$#" )); do
  case "$1" in
    -i) INPUT="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    -m) MODE="$2"; shift 2 ;;
    --report) REPORT=1; shift ;;
    --strict) STRICT=1; shift ;;
    -v|--verbose) VERBOSE=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --skip-complete) SKIP_COMPLETE=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

[[ -n "${INPUT}" && -n "${OUTDIR}" && -n "${MODE}" ]] || { usage; exit 1; }
[[ -f "${INPUT}" ]] || { echo "ERROR: input file not found: ${INPUT}" >&2; exit 1; }
[[ "${MODE}" == "symlink" || "${MODE}" == "concat" ]] || { echo "ERROR: -m must be 'symlink' or 'concat'" >&2; exit 1; }

log "Starting: mode=${MODE}, report=${REPORT}, strict=${STRICT}, verbose=${VERBOSE}, dry-run=${DRYRUN}, skip-complete=${SKIP_COMPLETE}, force=${FORCE}"
log "Input list: ${INPUT}"
log "Output dir: ${OUTDIR}"
doit "mkdir -p \"${OUTDIR}\""
TMPDIR="${OUTDIR}/.tmp_work"
doit "mkdir -p \"${TMPDIR}\""
TMP_R2="${TMPDIR}/autoinfer_r2.list"
[[ "${DRYRUN}" -eq 0 ]] && : > "${TMP_R2}" || true

# ---------- Data structures ----------
declare -A R1_LIST R2_LIST R1_LANES R2_LANES
declare -A SEEN_INPUT       # de-dupe the raw input lines
declare -A SEEN_TUBE_FILE   # de-dupe per tube+path (prevents double merge)
total_lines=0; valid_files=0; skipped_files=0; inferred_pairs=0

# ---------- Helpers ----------
parse_one_file() {
  local f="$1"
  [[ -n "${f}" ]] || return 0
  if [[ ! -f "${f}" ]]; then
    echo "WARN: missing file: ${f}" >&2; skipped_files=$((skipped_files+1)); return 0
  fi

  # Canonicalize for robust duplicate detection (handles symlinks/relative paths).
  local ap
  ap="$(readlink -f -- "$f" 2>/dev/null || realpath -- "$f" 2>/dev/null || echo "$f")"
  local fn; fn="$(basename -- "$ap")"

  # dna_tube_id: support _d###_ and _pos#_
  local tube=""
  if [[ "${fn}" =~ _d([0-9]+)_ ]]; then
    tube="d${BASH_REMATCH[1]}"
  elif [[ "${fn}" =~ _pos([0-9]+)_ ]]; then
    tube="pos${BASH_REMATCH[1]}"
  else
    echo "WARN: could not parse dna_tube_id (expected _d###_ or _pos#_) in: ${fn}; skipping" >&2
    skipped_files=$((skipped_files+1)); return 0
  fi

  # read token: _R1.fastq.gz or _R2.fastq.gz
  local readnum=""
  if [[ "${fn}" =~ _R([12])\.fastq\.gz$ ]]; then
    readnum="${BASH_REMATCH[1]}"
  else
    echo "WARN: could not parse read (R1/R2) in: ${fn}; skipping" >&2
    skipped_files=$((skipped_files+1)); return 0
  fi

  # lane token (optional): _L###_
  local lane="NA"
  if [[ "${fn}" =~ _L([0-9]{3})_ ]]; then
    lane="${BASH_REMATCH[1]}"
  fi

  # ---- NEW: per-tube de-duplication ----
  local key="${tube}|${ap}"
  if [[ -n "${SEEN_TUBE_FILE[$key]:-}" ]]; then
    vlog "Dedup: already queued for tube=${tube}: ${ap}"
    return 0
  fi
  SEEN_TUBE_FILE["$key"]=1

  if [[ "${readnum}" == "1" ]]; then
    R1_LIST["$tube"]+="${ap}"$'\n'
    R1_LANES["$tube"]+="${lane}"$'\n'
  else
    R2_LIST["$tube"]+="${ap}"$'\n'
    R2_LANES["$tube"]+="${lane}"$'\n'
  fi
  valid_files=$((valid_files+1))
  vlog "Parsed: tube=${tube} read=R${readnum} lane=${lane} file=${fn}"
}

sort_fastqs() {
  awk '
    {
      fname=$0
      key=fname
      match(fname, /_L([0-9]{3})_/, m)
      if (m[1] != "") { key = sprintf("%03d", m[1]) "_" fname } else { key = "999_" fname }
      print key "\t" fname
    }
  ' | sort -t$'\t' -k1,1V | cut -f2-
}

make_manifest() {
  # args: newline-separated file list
  awk 'NF' | while read -r p; do
    [[ -f "$p" ]] || continue
    ap="$(readlink -f "$p" 2>/dev/null || realpath "$p" 2>/dev/null || echo "$p")"
    sz="$(stat -c%s "$p" 2>/dev/null || stat -f%z "$p")"
    mt="$(date -r "$p" +%s 2>/dev/null || stat -f%m "$p")"
    printf "%s\t%s\t%s\n" "$ap" "$sz" "$mt"
  done | sort -u
}

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then
    md5 | awk '{print $4}'
  else
    cksum | awk '{print $1}'
  fi
}

newest_mtime() {
  awk 'NF' | while read -r p; do
    [[ -f "$p" ]] || continue
    date -r "$p" +%s 2>/dev/null || stat -f%m "$p"
  done | sort -n | tail -n1
}

outputs_up_to_date() {
  local tdir="$1"; shift
  local newest_in="$1"; shift
  local ok=1
  if [[ "${MODE}" == "concat" ]]; then
    for o in "${tdir}/${tube}_R1.fastq.gz" "${tdir}/${tube}_R2.fastq.gz"; do
      if [[ ! -f "$o" ]]; then ok=0; break; fi
      out_mt=$(date -r "$o" +%s 2>/dev/null || stat -f%m "$o")
      if [[ "$out_mt" -lt "$newest_in" ]]; then ok=0; break; fi
      if [[ ! -s "$o" ]]; then ok=0; break; fi
    done
  else
    :
  fi
  echo "$ok"
}

# ---------- Pass 1: read list & auto-infer R2 ----------
# NEW: de-dupe the raw input list up-front
while IFS= read -r f_raw; do
  total_lines=$((total_lines+1))
  [[ -n "${f_raw}" ]] || continue
  # normalize to absolute as the key for dedup
  f_abs="$(readlink -f -- "$f_raw" 2>/dev/null || realpath -- "$f_raw" 2>/dev/null || echo "$f_raw")"
  if [[ -n "${SEEN_INPUT[$f_abs]:-}" ]]; then
    vlog "Dedup(INPUT): ${f_abs}"
    continue
  fi
  SEEN_INPUT["$f_abs"]=1

  # If R1, enqueue R2 only if not already present in input set
  if [[ "${f_abs}" =~ _R1\.fastq\.gz$ ]]; then
    r2="${f_abs/_R1.fastq.gz/_R2.fastq.gz}"
    if [[ -f "${r2}" && -z "${SEEN_INPUT[$r2]:-}" ]]; then
      inferred_pairs=$((inferred_pairs+1))
      [[ "${DRYRUN}" -eq 0 ]] && echo "${r2}" >> "${TMP_R2}" || echo "DRY-RUN: would auto-enqueue R2: ${r2}"
      SEEN_INPUT["$r2"]=1
    elif [[ ! -f "${r2}" ]]; then
      echo "WARN: no R2 found for R1: ${f_abs}" >&2
    fi
  fi

  parse_one_file "${f_abs}"
done < <(sed 's/\r$//' "${INPUT}" | awk 'NF' | sort -u)

log "Pass 1 complete: ${total_lines} lines read | ${valid_files} files parsed | ${inferred_pairs} R2 mates enqueued | ${skipped_files} skipped."

# Pass 2: ingest auto-inferred R2s
if [[ -s "${TMP_R2}" ]]; then
  vlog "Ingesting auto-inferred R2 list from ${TMP_R2}"
  # TMP_R2 is already deduped via SEEN_INPUT above, but de-dupe again for safety
  while IFS= read -r f; do parse_one_file "${f}"; done < <(sort -u "${TMP_R2}")
fi

# ---------- Reporting setup ----------
REPORT_DIR="${OUTDIR}/report"
if [[ "${REPORT}" -eq 1 ]]; then
  doit "mkdir -p \"${REPORT_DIR}\""
  TSV="${REPORT_DIR}/fastq_lane_summary_by_dna_tube_id.tsv"
  CSV="${REPORT_DIR}/fastq_lane_summary_by_dna_tube_id.csv"
  [[ "${DRYRUN}" -eq 0 ]] && printf "dna_tube_id\ttotal_fastq_files\tn_R1\tn_R2\tn_unique_lanes_R1\tn_unique_lanes_R2\tlanes_R1\tlanes_R2\tstatus\n" > "${TSV}" || true
fi

# ---------- Process per tube ----------
STRICT_FAIL=0
tubes_processed=0
tubes_skipped_complete=0

for tube in $(printf "%s\n%s\n" "${!R1_LIST[@]}" "${!R2_LIST[@]}" | sort -u); do
  [[ -n "${tube}" ]] || continue
  tubes_processed=$((tubes_processed+1))
  tdir="${OUTDIR}/${tube}"
  doit "mkdir -p \"${tdir}\""

  # Collect & sort; NEW: final de-dup here too
  r1_files="$(printf "%s" "${R1_LIST[$tube]-}" | grep -v '^$' || true)"
  r2_files="$(printf "%s" "${R2_LIST[$tube]-}" | grep -v '^$' || true)"
  r1_sorted="$(printf "%s\n" ${r1_files} | grep -v '^$' | sort_fastqs | sort -u || true)"
  r2_sorted="$(printf "%s\n" ${r2_files} | grep -v '^$' | sort_fastqs | sort -u || true)"

  lanes_R1="$(printf "%s" "${R1_LANES[$tube]-}" | grep -v -e '^$' -e '^NA$' | sort -u | paste -sd, - || true)"
  lanes_R2="$(printf "%s" "${R2_LANES[$tube]-}" | grep -v -e '^$' -e '^NA$' | sort -u | paste -sd, - || true)"
  n_lanes_R1=0; [[ -n "${lanes_R1}" ]] && n_lanes_R1=$(tr ',' '\n' <<< "${lanes_R1}" | wc -l | tr -d ' ')
  n_lanes_R2=0; [[ -n "${lanes_R2}" ]] && n_lanes_R2=$(tr ',' '\n' <<< "${lanes_R2}" | wc -l | tr -d ' ')
  n_R1=0; [[ -n "${r1_sorted}" ]] && n_R1=$(wc -w <<< ${r1_sorted} | tr -d ' ')
  n_R2=0; [[ -n "${r2_sorted}" ]] && n_R2=$(wc -w <<< ${r2_sorted} | tr -d ' ')
  total=$(( n_R1 + n_R2 ))

  log "▶ [${tubes_processed}] ${tube}: files=${total} (R1=${n_R1}, R2=${n_R2}), lanes (R1=${n_lanes_R1}:[${lanes_R1:-NA}] | R2=${n_lanes_R2}:[${lanes_R2:-NA}])"

  # ----- Completion detection -----
  cur_manifest="$( (printf "%s\n" ${r1_sorted}; printf "%s\n" ${r2_sorted}) | make_manifest )"
  cur_hash="$(printf "%s" "${cur_manifest}" | md5_of)"
  newest_in="$( (printf "%s\n" ${r1_sorted}; printf "%s\n" ${r2_sorted}) | newest_mtime )"
  man_path="${tdir}/.inputs.manifest.txt"
  hash_path="${tdir}/.inputs.manifest.md5"

  up_to_date=0
  if [[ "${MODE}" == "concat" ]]; then
    up_to_date=$(outputs_up_to_date "${tdir}" "${newest_in:-0}")
  else
    up_to_date=1
  fi

  prev_hash=""
  if [[ -f "${hash_path}" ]]; then
    prev_hash="$(cat "${hash_path}" 2>/dev/null || true)"
  fi

  if [[ "${SKIP_COMPLETE}" -eq 1 && "${FORCE}" -eq 0 && "${up_to_date}" -eq 1 && -n "${prev_hash}" && "${prev_hash}" == "${cur_hash}" ]]; then
    log "   ✓ Skipping '${tube}' (already completed: outputs up-to-date & manifest unchanged)"
    tubes_skipped_complete=$((tubes_skipped_complete+1))
    status="OK_SKIPPED"
    if [[ "${REPORT}" -eq 1 && "${DRYRUN}" -eq 0 ]]; then
      printf "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" \
        "${tube}" "${total}" "${n_R1}" "${n_R2}" "${n_lanes_R1}" "${n_lanes_R2}" \
        "${lanes_R1:-NA}" "${lanes_R2:-NA}" "${status}" >> "${TSV}"
    fi
    continue
  fi

  # ----- Usual integrity checks -----
  status="OK"
  if [[ ${n_R1} -eq 0 || ${n_R2} -eq 0 ]]; then
    status="MISSING_R1_OR_R2"
    [[ "${STRICT}" -eq 1 ]] && STRICT_FAIL=1
  elif [[ ${n_lanes_R1} -ne ${n_lanes_R2} ]]; then
    status="LANE_COUNT_MISMATCH"
    [[ "${STRICT}" -eq 1 ]] && STRICT_FAIL=1
  fi

  # ----- Perform action -----
  if [[ "${MODE}" == "symlink" ]]; then
    if [[ -n "${r1_sorted}" ]]; then
      while IFS= read -r f; do [[ -n "$f" ]] && doit "ln -sf \"$(realpath "$f")\" \"${tdir}/\""; done <<< "${r1_sorted}"
    fi
    if [[ -n "${r2_sorted}" ]]; then
      while IFS= read -r f; do [[ -n "$f" ]] && doit "ln -sf \"$(realpath "$f")\" \"${tdir}/\""; done <<< "${r2_sorted}"
    fi
  else
    if [[ -n "${r1_sorted}" ]]; then
      out1="${tdir}/${tube}_R1.fastq.gz"
      log "  Concatenate R1 -> ${out1}"
      if [[ "${DRYRUN}" -eq 1 ]]; then
        echo "DRY-RUN: would cat R1 list into ${out1}"
      else
        r1_listfile="${TMPDIR}/${tube}_R1.list"
        printf "%s\n" ${r1_sorted} | sort -u > "${r1_listfile}"
        xargs -a "${r1_listfile}" -r cat > "${out1}"
      fi
    else
      echo "WARN: no R1 files for ${tube}" >&2
    fi
    if [[ -n "${r2_sorted}" ]]; then
      out2="${tdir}/${tube}_R2.fastq.gz"
      log "  Concatenate R2 -> ${out2}"
      if [[ "${DRYRUN}" -eq 1 ]]; then
        echo "DRY-RUN: would cat R2 list into ${out2}"
      else
        r2_listfile="${TMPDIR}/${tube}_R2.list"
        printf "%s\n" ${r2_sorted} | sort -u > "${r2_listfile}"
        xargs -a "${r2_listfile}" -r cat > "${out2}"
      fi
    else
      echo "WARN: no R2 files for ${tube}" >&2
    fi
  fi

  # ----- Write/update manifest + hash (if not dry-run) -----
  if [[ "${DRYRUN}" -eq 0 ]]; then
    printf "%s\n" "${cur_manifest}" > "${man_path}"
    printf "%s\n" "${cur_hash}" > "${hash_path}"
  fi

  # ----- Report row -----
  if [[ "${REPORT}" -eq 1 && "${DRYRUN}" -eq 0 ]]; then
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" \
      "${tube}" "${total}" "${n_R1}" "${n_R2}" "${n_lanes_R1}" "${n_lanes_R2}" \
      "${lanes_R1:-NA}" "${lanes_R2:-NA}" "${status}" >> "${TSV}"
  fi
done

if [[ "${REPORT}" -eq 1 && "${DRYRUN}" -eq 0 ]]; then
  awk -F'\t' 'BEGIN{OFS=","} {gsub(/\t/,","); print}' "${TSV}" > "${CSV}"
  log "Report written: ${TSV}"
  log "Report written: ${CSV}"
fi

if [[ "${STRICT}" -eq 1 && "${STRICT_FAIL}" -ne 0 ]]; then
  echo "ERROR: strict checks failed (see status in report/log)." >&2
  exit 2
fi

log "Done. Tubes processed: ${tubes_processed}. Skipped (complete): ${tubes_skipped_complete}. Output under: ${OUTDIR}"
