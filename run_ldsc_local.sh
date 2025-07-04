#!/usr/bin/env bash
set -euo pipefail
set -x

## 0) Sanity‐check tools
which munge_sumstats.py
munge_sumstats.py --help | head -n2
which ldsc.py
ldsc.py --help | head -n2

## 1) Paths
WORKDIR="$HOMEpath/cardiovascularNewAnalysis/new_LDSC_analysis"
cd "$WORKDIR"

## 2) HapMap3 SNP‐list
MERGE_ALLELES="$WORKDIR/w_hm3.snplist"
[[ -f "$MERGE_ALLELES" ]] || { echo "ERROR: missing $MERGE_ALLELES"; exit 1; }

## 3) LD-score prefixes
#   EUR: the standard eur_w_ld_chr/1.l2.ldscore.gz … 22.l2.ldscore.gz  
#   AFR: your ALL_ensemble…_chr1.l2.ldscore.gz … chr22.l2.ldscore.gz
declare -A LD_PREFIX=(
  [EUR]="$WORKDIR/ldsc_data/eur_w_ld_chr/"
  [AFR]="$WORKDIR/ldsc_data/AFR/ALL_ensemble_1000G_hg38_AFR_chr@"
)

# sanity checks
[[ -d "${LD_PREFIX[EUR]}" ]] \
  || { echo "ERROR: missing EUR panel at ${LD_PREFIX[EUR]}"; exit 1; }

AFR_CHR1="${WORKDIR}/ldsc_data/AFR/ALL_ensemble_1000G_hg38_AFR_chr1.l2.ldscore.gz"
[[ -f "$AFR_CHR1" ]] \
  || { echo "ERROR: missing AFR chr1 file at $AFR_CHR1"; exit 1; }

## 4) Traits & sample sizes
TRAITS=(SystolicBP DiastolicBP PulseRate MaxHeartRate)
declare -A N_SAMPLE=( [AFR]=6551 [EUR]=396670 )
ANCESTRIES=(AFR EUR)

## 5) Munge (if you haven’t already)
for anc in "${ANCESTRIES[@]}"; do
 for trait in "${TRAITS[@]}"; do
   munge_sumstats.py \
     --sumstats "ldsc_${trait}_${anc}.sumstats.txt" \
     --N "${N_SAMPLE[$anc]}" \
     --out "munged_${trait}_${anc}" \
     --merge-alleles "$MERGE_ALLELES"
 done
done

## 6) LDSC h² & rg
for anc in "${ANCESTRIES[@]}"; do
  echo
  echo "=== HERITABILITY for $anc ==="
  for trait in "${TRAITS[@]}"; do
    SUMSTATS="munged_${trait}_${anc}.sumstats.gz"
    OUT="${trait}_${anc}_h2"
    ldsc.py \
      --h2 "$SUMSTATS" \
      --ref-ld-chr "${LD_PREFIX[$anc]}" \
      --w-ld-chr   "${LD_PREFIX[$anc]}" \
      --out "$OUT"
  done

  echo
  echo "=== GENETIC CORRELATION for $anc ==="
  RG_LIST=""
  for trait in "${TRAITS[@]}"; do
	RG_LIST+="munged_${trait}_${anc}.sumstats.gz,"
  done
  RG_LIST="${RG_LIST%,}"
	
  echo "→ ldsc.py --rg $RG_LIST → ${anc}_rg"
  ldsc.py \
	--rg "$RG_LIST" \
	--ref-ld-chr "${LD_PREFIX[$anc]}" \
	--w-ld-chr   "${LD_PREFIX[$anc]}" \
	--out "${anc}_rg"
done

echo "✅ LDSC pipeline completed locally!"
