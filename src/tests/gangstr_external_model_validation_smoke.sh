#!/usr/bin/env bash
set -euo pipefail

gangstr_bin="$1"
repo_root="$2"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

regions="$tmpdir/CACNA1A_local.tsv"
printf '19\t5000\t5039\t3\tCTG\n' > "$regions"

common_args=(
  --bam "$repo_root/unused/tests/54_nc_12.sorted.bam"
  --ref "$repo_root/unused/tests/CACNA1A_5k_region.fa"
  --regions "$regions"
  --coverage 80
  --insertmean 500
  --insertsdev 50
  --readlength 100
  --stutter-mode external
)

cat > "$tmpdir/good.models.tsv" <<'EOF'
chrom	start	end	period	motif	in_geom	in_up	in_down	out_geom	out_up	out_down
19	5000	5039	3	ctg	0.982278	0.030949	0.113803	0.953461	0.00915561	0.00915561
EOF

"$gangstr_bin" "${common_args[@]}" \
  --stutter-model-tsv "$tmpdir/good.models.tsv" \
  --out "$tmpdir/good" \
  --verbose > "$tmpdir/good.stdout" 2> "$tmpdir/good.stderr"
grep -q "Loaded 1 external stutter models" "$tmpdir/good.stderr"
test -s "$tmpdir/good.vcf"

cat > "$tmpdir/bad_numeric.models.tsv" <<'EOF'
chrom	start	end	period	motif	in_geom	in_up	in_down	out_geom	out_up	out_down
19	5000	5039	3	ctg	0.982278	0.030949junk	0.113803	0.953461	0.00915561	0.00915561
EOF

if "$gangstr_bin" "${common_args[@]}" \
    --stutter-model-tsv "$tmpdir/bad_numeric.models.tsv" \
    --out "$tmpdir/bad_numeric" \
    > "$tmpdir/bad_numeric.stdout" 2> "$tmpdir/bad_numeric.stderr"; then
  echo "Expected bad numeric external model to fail" >&2
  exit 1
fi
grep -q "not a valid finite number" "$tmpdir/bad_numeric.stderr"

cat > "$tmpdir/missing_column.models.tsv" <<'EOF'
chrom	start	end	period	motif	in_geom	in_down	out_geom	out_up	out_down	junk
19	5000	5039	3	ctg	0.982278	0.113803	0.953461	0.00915561	0.00915561	0.01
EOF

if "$gangstr_bin" "${common_args[@]}" \
    --stutter-model-tsv "$tmpdir/missing_column.models.tsv" \
    --out "$tmpdir/missing_column" \
    > "$tmpdir/missing_column.stdout" 2> "$tmpdir/missing_column.stderr"; then
  echo "Expected missing-column external model to fail" >&2
  exit 1
fi
grep -q "missing required column: in_up" "$tmpdir/missing_column.stderr"

cat > "$tmpdir/up_gt_down.models.tsv" <<'EOF'
chrom	start	end	period	motif	in_geom	in_up	in_down	out_geom	out_up	out_down
19	5000	5039	3	ctg	0.982278	0.20	0.05	0.953461	0.00915561	0.00915561
EOF

if "$gangstr_bin" "${common_args[@]}" \
    --stutter-model-tsv "$tmpdir/up_gt_down.models.tsv" \
    --out "$tmpdir/up_gt_down" \
    > "$tmpdir/up_gt_down.stdout" 2> "$tmpdir/up_gt_down.stderr"; then
  echo "Expected in_up > in_down external model to fail" >&2
  exit 1
fi
grep -q "P_UP must not exceed P_DOWN" "$tmpdir/up_gt_down.stderr"
