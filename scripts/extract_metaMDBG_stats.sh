#!/usr/bin/env bash
# ------------------------------------------------------------
# Script: extract_metaMDBG_stats.sh
# Purpose: Extract assembly statistics from metaMDBG.log files
# ------------------------------------------------------------

set -euo pipefail

# Base directory containing assemblies
base="/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/lr_results/assemblies"

# Output summary file
out="metaMDBG_assembly_summary.tsv"

# Write header
echo -e "sample\tassembly_length\tcontigs_n50\tnb_contigs\tnb_contigs_gt1mb\tnb_circular_contigs_gt1mb" > "$out"

# Search:
# 1) Top-level assemblies (coassemblies, hybrid, ONT)
# 2) HiFi individual assemblies
for d in "$base"/* "$base"/hifi_individual/*; do
  [[ -d "$d" ]] || continue
  [[ -f "$d/contigs.fasta.gz" ]] || continue
  [[ -f "$d/metaMDBG.log" ]] || continue

  sample="$(basename "$d")"

  # Parse the LAST occurrence of each stat from metaMDBG.log
  awk -v sample="$sample" '
    BEGIN {
      al=n50=nc=gt1=circ=""
    }

    /^\s*Assembly length:/ {
      al=$0
      sub(/.*Assembly length:\s*/, "", al)
      gsub(/[^0-9]/, "", al)
    }

    /^\s*Contigs N50:/ {
      n50=$0
      sub(/.*Contigs N50:\s*/, "", n50)
      gsub(/[^0-9]/, "", n50)
    }

    /^\s*Nb contigs:/ {
      nc=$0
      sub(/.*Nb contigs:\s*/, "", nc)
      gsub(/[^0-9]/, "", nc)
    }

    /^\s*Nb Contigs \(>1Mb\):/ {
      gt1=$0
      sub(/.*Nb Contigs \(>1Mb\):\s*/, "", gt1)
      gsub(/[^0-9]/, "", gt1)
    }

    /^\s*Nb circular contigs \(>1Mb\):/ {
      circ=$0
      sub(/.*Nb circular contigs \(>1Mb\):\s*/, "", circ)
      gsub(/[^0-9]/, "", circ)
    }

    END {
      if (al   == "") al   = "NA"
      if (n50  == "") n50  = "NA"
      if (nc   == "") nc   = "NA"
      if (gt1  == "") gt1  = "NA"
      if (circ == "") circ = "NA"

      printf "%s\t%s\t%s\t%s\t%s\t%s\n", sample, al, n50, nc, gt1, circ
    }
  ' "$d/metaMDBG.log" >> "$out"

done

echo "Wrote: $out"
