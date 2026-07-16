#!/bin/bash
## ============================================================
## submit_all_caqtl.sh
## Submits all caQTL jobs in batches respecting SLURM MaxArraySize
##
## Usage:
##   bash submit_all_caqtl.sh gemma
##   bash submit_all_caqtl.sh gemmamod
##
## Automatically calculates number of batches needed and submits
## all of them, printing the job ID for each.
## ============================================================

PEAK_LIST="caqtl_gene_inputs/peak_list.txt"
ARRAY_SCRIPT="submit_gemma_array_03.sh"
MAX_ARRAY=5000    # stay safely under the 20000 limit
CONCURRENT=375   # max jobs running at once per batch

if [ ! -f "$PEAK_LIST" ]; then
  echo "ERROR: Peak list not found at $PEAK_LIST"
  exit 1
fi

TOTAL=$(wc -l < "$PEAK_LIST")
echo "Total peaks     : $TOTAL"
echo "Batch size      : $MAX_ARRAY"
echo "Max concurrent  : $CONCURRENT"
echo ""

# Calculate number of batches needed
N_BATCHES=$(( (TOTAL + MAX_ARRAY - 1) / MAX_ARRAY ))
echo "Batches needed  : $N_BATCHES"
echo "-----------------------------------"

PREV_JOB_ID=""

for BATCH in $(seq 1 $N_BATCHES); do

  OFFSET=$(( (BATCH - 1) * MAX_ARRAY ))
  BATCH_START=1
  BATCH_END=$(( BATCH * MAX_ARRAY ))

  # Last batch may be smaller
  if [ $BATCH_END -gt $TOTAL ]; then
    BATCH_END=$TOTAL
  fi

  BATCH_SIZE=$(( BATCH_END - OFFSET ))

  echo "Submitting batch $BATCH / $N_BATCHES"
  echo "  Peaks ${OFFSET+1} – ${BATCH_END} (${BATCH_SIZE} peaks, offset=${OFFSET})"
  
  if [ -z "$PREV_JOB_ID" ]; then
 
    # First batch: no dependency, submit immediately
    JOB_ID=$(sbatch \
      --account=pro00117530 \
      --array=1-${BATCH_SIZE}%${CONCURRENT} \
      --export=ALL,OFFSET=${OFFSET} \
      --job-name=caqtl_b${BATCH} \
      --output=logs/caqtl_b${BATCH}_%A_%a.out \
      --error=logs/caqtl_b${BATCH}_%A_%a.err \
      --cpus-per-task=1 \
      --mem=1G \
      --time=04:00:00 \
      "$ARRAY_SCRIPT" \
      | awk '{print $NF}')
 
    echo "  Submitted (no dependency) — job ID: $JOB_ID"
 
  else
 
    # All subsequent batches: hard dependency on previous batch
    JOB_ID=$(sbatch \
      --account=pro00117530 \
      --dependency=afterok:${PREV_JOB_ID} \
      --array=1-${BATCH_SIZE}%${CONCURRENT} \
      --export=ALL,OFFSET=${OFFSET} \
      --job-name=caqtl_b${BATCH} \
      --output=logs/caqtl_b${BATCH}_%A_%a.out \
      --error=logs/caqtl_b${BATCH}_%A_%a.err \
      --cpus-per-task=1 \
      --mem=1G \
      --time=04:00:00 \
      "$ARRAY_SCRIPT" \
      | awk '{print $NF}')
 
    echo "  Submitted (depends on job $PREV_JOB_ID) — job ID: $JOB_ID"
 
  fi
 
  # Validate that sbatch actually returned a job ID
  if ! [[ "$JOB_ID" =~ ^[0-9]+$ ]]; then
    echo "ERROR: sbatch did not return a valid job ID for batch $BATCH."
    echo "  Got: '$JOB_ID'"
    echo "  Stopping submission to avoid unthrottled queue flooding."
    exit 1
  fi
 
  PREV_JOB_ID=$JOB_ID
  echo ""
  sleep 2
 
done
 
echo "All $N_BATCHES batches queued sequentially."
echo "They will run ONE AT A TIME — each waits for the previous to finish."
echo ""
echo "Monitor queue:"
echo "  squeue -u $USER -o '%.18i %.30j %.8T %.10M' | grep caqtl"
echo ""
echo "Check output progress:"
echo "  ls caqtl_gene_outputs/*.rds 2>/dev/null    | grep -v ERROR | wc -l"
echo "  ls caqtl_gemmamod_outputs/*.rds 2>/dev/null | grep -v ERROR | wc -l"
echo ""
echo "To cancel everything:"
echo "  scancel \$(squeue -u $USER -h -o '%i' | tr '\n' ' ')"
