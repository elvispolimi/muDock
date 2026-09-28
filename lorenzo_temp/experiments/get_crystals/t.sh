ids=(1fkb 1hii 2ya6 3udd 4few 5cst 5uez 5wuk)
crystal_scores=(-15.797226 -15.021526 -10.792323 -18.935410 244.454712 -3.826046 -9.453270 -10.575790)

for i in "${!ids[@]}"; do
    id="${ids[$i]}"
    crystal_score="${crystal_scores[$i]}"

    echo "$id $crystal_score"
done

for dir in ${DATA_DIR}/*/; do