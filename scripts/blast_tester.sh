#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
BLAST remote tester / runner.

Usage:
  scripts/blast_tester.sh [options]

Options:
  --query <fasta>           Query FASTA. If omitted, uses a built-in test query.
  --out-prefix <path>       Output prefix (default: results/blast_test)
  --db <name>               BLAST database (default: nt)
  --task <task>             blastn task (default: megablast)
  --evalue <val>            E-value cutoff (default: 1e-20)
  --max-target-seqs <n>     Max target seqs (default: 20)
  --timeout-sec <n>         Timeout in seconds (default: 900)
  --no-remote               Do not use -remote (requires local BLAST DB)
  -h, --help                Show this help

Examples:
  scripts/blast_tester.sh
  scripts/blast_tester.sh --query results/LFFHG8_1_pcr1_sub5/LFFHG8_1_pcr1_sub5.consensus.fasta \
    --out-prefix results/LFFHG8_1_pcr1_sub5/contaminant_blast
USAGE
}

QUERY=""
OUT_PREFIX="results/blast_test"
DB="nt"
TASK="megablast"
EVALUE="1e-20"
MAX_TARGET="20"
TIMEOUT_SEC="900"
USE_REMOTE="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --query)
      QUERY="$2"; shift 2 ;;
    --out-prefix)
      OUT_PREFIX="$2"; shift 2 ;;
    --db)
      DB="$2"; shift 2 ;;
    --task)
      TASK="$2"; shift 2 ;;
    --evalue)
      EVALUE="$2"; shift 2 ;;
    --max-target-seqs)
      MAX_TARGET="$2"; shift 2 ;;
    --timeout-sec)
      TIMEOUT_SEC="$2"; shift 2 ;;
    --no-remote)
      USE_REMOTE="0"; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1 ;;
  esac
done

if ! command -v blastn >/dev/null 2>&1; then
  echo "blastn not found in PATH" >&2
  exit 1
fi

TMP_QUERY=""
if [[ -z "$QUERY" ]]; then
  TMP_QUERY="$(mktemp -t blast_test_query.XXXXXX.fa)"
  cat > "$TMP_QUERY" <<'QEOF'
>test_query_16S_like
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAACTGAGA
QEOF
  QUERY="$TMP_QUERY"
  echo "No query provided; using built-in test query: $QUERY"
fi

cleanup() {
  if [[ -n "$TMP_QUERY" && -f "$TMP_QUERY" ]]; then
    rm -f "$TMP_QUERY"
  fi
}
trap cleanup EXIT

if [[ ! -f "$QUERY" ]]; then
  echo "Query file not found: $QUERY" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"
OUT_TSV="${OUT_PREFIX}.tsv"
OUT_LOG="${OUT_PREFIX}.log"

OUTFMT="6 qseqid sacc pident length qcovs qstart qend sstart send evalue bitscore sscinames stitle"

CMD=(blastn)
if [[ "$USE_REMOTE" == "1" ]]; then
  CMD+=(-remote)
fi
CMD+=(
  -query "$QUERY"
  -db "$DB"
  -task "$TASK"
  -evalue "$EVALUE"
  -max_target_seqs "$MAX_TARGET"
  -outfmt "$OUTFMT"
  -out "$OUT_TSV"
)

{
  echo "timestamp: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "blastn_version: $(blastn -version | head -n 1)"
  echo "query: $QUERY"
  echo "out_tsv: $OUT_TSV"
  echo "db: $DB"
  echo "task: $TASK"
  echo "remote: $USE_REMOTE"
  echo "evalue: $EVALUE"
  echo "max_target_seqs: $MAX_TARGET"
  echo "timeout_sec: $TIMEOUT_SEC"
  echo "command: ${CMD[*]}"
} > "$OUT_LOG"

set +e
if timeout "$TIMEOUT_SEC" "${CMD[@]}" >> "$OUT_LOG" 2>&1; then
  rc=0
else
  rc=$?
fi
set -e

lines=0
if [[ -f "$OUT_TSV" ]]; then
  lines=$(wc -l < "$OUT_TSV")
fi

echo "BLAST exit code: $rc"
echo "Output TSV: $OUT_TSV"
echo "Output LOG: $OUT_LOG"
echo "Hit lines: $lines"

if [[ "$rc" -eq 124 ]]; then
  echo "Timed out before completion. Increase --timeout-sec or check network access to NCBI BLAST." >&2
  exit 124
fi

if [[ "$rc" -ne 0 ]]; then
  echo "BLAST failed. See log: $OUT_LOG" >&2
  exit "$rc"
fi

if [[ "$lines" -eq 0 ]]; then
  echo "No hits returned (or output empty)." >&2
  exit 2
fi

echo "Top hits:"
sed -n '1,10p' "$OUT_TSV"
