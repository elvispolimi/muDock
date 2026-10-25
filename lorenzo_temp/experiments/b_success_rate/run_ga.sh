#!/usr/bin/env bash
IGNORE_TOKEN="Experiment"

MAX_JOBS=12

BUILD=build
DATA_DIR="./data/coreset_CASF_2016"
EXP_DIR="./lorenzo_temp/experiments"
TEST_DIR="${EXP_DIR}/test"

ids=(1a30 1bzc 1c5z 1e66 1eby 1g2k 1gpk 1gpn 1h22 1h23 1k1i 1lpg 1mq6 1nc1 1nc3 1nvq 1o0h 1o3f 1o5b 1owh 1p1n 1p1q 1ps3 1pxn 1q8t 1q8u 1qf1 1qkt 1r5y 1s38 1sqa 1syi 1u1b 1uto 1vso 1w4o 1y6r 1yc1 1ydr 1ydt 1z6e 1z95 1z9g 2al5 2br1 2brb 2c3i 2cbv 2cet 2fvd 2fxs 2hb1 2iwx 2j78 2j7h 2p15 2p4y 2pog 2qbp 2qbq 2qbr 2qe4 2qnq 2v00 2v7a 2vkm 2vvn 2vw5 2w4x 2w66 2wbg 2wca 2weg 2wer 2wn9 2wnc 2wtv 2wvt 2x00 2xb8 2xbv 2xdl 2xii 2xj7 2xnb 2xys 2y5h 2yfe 2yge 2yki 2ymd 2zb1 2zcq 2zcr 2zy1 3acw 3ag9 3arp 3arq 3aru 3arv 3ary 3b1m 3b27 3b5r 3b65 3b68 3bgz 3cj4 3coy 3coz 3d4z 3d6q 3dd0 3dx1 3dx2 3dxg 3e5a 3e92 3e93 3ebp 3ejr 3f3a 3f3c 3f3d 3f3e 3fcq 3fur 3fv1 3fv2 3g0w 3g2n 3g2z 3g31 3gbb 3gc5 3ge7 3gnw 3gy4 3ivg 3jvr 3jvs 3jya 3k5v 3kgp 3kr8 3kwa 3l7b 3mss 3myg 3n76 3n7a 3n86 3nq9 3nw9 3o9i 3oe4 3oe5 3ozs 3ozt 3p5o 3prs 3pww 3pxf 3pyy 3qgy 3qqs 3r88 3rlr 3rr4 3rsx 3ryj 3syr 3twp 3u5j 3u8k 3u8n 3u9q 3udh 3ueu 3uev 3uew 3uex 3ui7 3uo4 3up2 3uri 3uuo 3wtj 3wz8 3zdg 4abg 4agn 4agp 4agq 4bkt 4ciw 4cr9 4cra 4crc 4ddh 4ddk 4de1 4de2 4de3 4djv 4dld 4dli 4e5w 4e6q 4ea2 4eky 4eo8 4eor 4f09 4f2w 4f3c 4f9w 4gfm 4gid 4gkm 4hge 4ih5 4ih7 4ivb 4ivc 4ivd 4j21 4j28 4j3l 4jfs 4jia 4jsz 4k18 4k77 4kzq 4kzu 4llx 4lzs 4m0y 4m0z 4mgd 4mme 4ogj 4owm 4pcs 4qac 4qd6 4rfm 4tmn 4twp 4ty7 4u4s 4w9c 4w9h 4w9i 4w9l 4wiv 4x6p 5a7b 5aba 5c28 5c2h 5dwr 5tmn)
crystal_scores=(-4.92767 -8.844632 -0.851277 -6.59583 -8.596169 -7.227155 -1.729623 0.339111 -7.53787 288.158691 -3.585232 -1.29267 -3.988441 -2.490548 -1.866436 -3.18403 4.153614 41.424843 -2.887584 -2.705627 -6.042484 -4.387009 -3.324558 0.72043 -2.390844 -2.615715 -10.704685 -3.122935 -1.600831 -3.51778 0.116504 -8.097101 8.905769 -2.418868 -5.97896 0.571738 -2.329404 -3.772542 -3.484977 -3.822094 -3.577033 -4.277293 -3.668527 -7.630584 -2.72521 -3.925623 -1.101734 -1.399765 -0.67478 -2.493005 -1.407803 -3.928283 -3.609652 -2.097933 0.456109 -4.373397 -4.431096 -4.337811 -1.386554 -5.157811 -4.475329 -3.012466 -2.983602 -2.363991 -4.851544 -7.217653 -1.182074 -3.711936 -1.177679 -2.762384 -2.19136 -2.046038 -2.241997 5.053053 -2.716562 -2.064045 -5.167174 -2.177783 -2.220733 -5.933784 -4.623016 -4.274053 -2.897148 -2.412557 -1.919087 2.693094 -1.297716 9.419469 3.386631 -9.488885 -3.527147 -6.22124 -3.441498 -6.819621 -3.807363 -2.794099 5.16925 -8.735242 78.44294 -2.352705 56044.804688 2.897366 -4.314973 -2.65517 -4.727698 -6.134722 -2.778485 -5.351007 -2.099269 -4.620814 -3.35963 -3.510527 -0.165224 21.476011 -3.659555 -5.2434 -1.734024 31872.945312 -7.909467 -8.255076 -1.080359 -7.928732 -3.969706 -5.518228 -4.488379 -5.911496 -4.593515 -7.530247 -6.779161 -8.291726 -2.11481 -0.515949 -1.790968 -0.285482 -7.288947 -3.483639 -1.92139 -4.973276 0.067972 -0.101351 -3.832939 -6.613667 -2.213093 -2.059274 -1.045209 -4.869537 74060.648438 1.151293 -2.22304 -3.169144 -4.182385 -2.825772 -5.44331 -0.63603 -8.299557 -6.878663 -9.97984 -4.006662 -0.747438 -1.960122 -5.091401 -6.458052 -5.720702 0.127209 -3.734506 -2.38641 -2.376544 -1.099299 -5.002181 -2.591565 -3.309969 -4.416487 -1.640649 24.532885 -4.21227 -2.050527 -3.763513 -3.705533 -1.451048 -2.368228 -2.921494 -3.228291 -3.934497 -3.50194 -3.850298 -4.279899 -9.374928 -3.639116 -1.827208 -6.445173 -1.133021 -0.482512 3.680776 0.161638 -0.544027 0.974669 0.460841 -2.532456 1.062327 -8.798751 -4.653146 -2.268387 -4.927009 -4.036175 -5.847238 -0.776836 -7.180317 -4.013467 -7.283812 -2.633337 1314.692627 -0.925193 -5.736087 -3.717469 -3.118294 -1.238 -2.818328 -5.927937 -3.722718 -4.41068 -3.134285 -4.370174 146.485519 -3.538471 -5.919197 -5.618873 -6.703081 -6.526631 -2.562396 -7.419788 -1.080559 0.227623 0.070017 -3.530325 -2.843251 -6.22261 -8.078053 -1.868324 -3.760489 0.571722 -5.502382 7.306289 -4.990584 -1.042561 -1.19959 -2.203516 -1.546753 1.998976 1.362901 3.588154 111.937057 -8.901407 0.335885 17.107246 36.803482 90.848541 136.98555 -5.579881 -8.5956 -1.606827 4.755371 -3.138379 -9.178081 -6.155669 -11.531065)

NUM_SEEDS=10
seeds=($(seq 0 $((NUM_SEEDS - 1))))

POPULATION=100
GENERATIONS=1000000 # high value => not bounded
BUDGET=30 #seconds

TOTAL_RUNS=$((${#ids[@]} * ${#seeds[@]}))
COMPLETED_RUNS=0

rm ${TEST_DIR}/*
echo "Docking..."
for SEED in "${seeds[@]}"; do
    for i in "${!ids[@]}"; do
        PDBID="${ids[$i]}"
        CRYSTAL_SCORE="${crystal_scores[$i]}"

        POCKET="${DATA_DIR}/${PDBID}/${PDBID}_pocket.pdb"
        LIGAND="${DATA_DIR}/${PDBID}/${PDBID}_ligand.adtmol2"        
        
        OUT_PATH="${TEST_DIR}/${PDBID}_${LSRATE}_${LSIT}_${SEED}.txt"

        ./builds/"$BUILD"/application/muDock \
            --protein "$POCKET" \
            --ligand "$LIGAND" \
            --seed "$SEED" \
            --search genetic \
            --population "$POPULATION" \
            --generations "$GENERATIONS" \
            --crystal_score "$CRYSTAL_SCORE" \
            --placement preserve \
            --max-search-box 5 \
            --max_docking_time "$BUDGET" \
            > "$OUT_PATH" 2>&1 &

        if [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; then
            wait -n
            COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
            printf "\rProgress: [%d/%d] %3d%%" \
                "$COMPLETED_RUNS" \
                "$TOTAL_RUNS" \
                "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
        fi

    done
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
    wait -n
    COMPLETED_RUNS=$((COMPLETED_RUNS + 1))
    printf "\rProgress: [%d/%d] %3d%%" \
        "$COMPLETED_RUNS" \
        "$TOTAL_RUNS" \
        "$((COMPLETED_RUNS * 100 / TOTAL_RUNS))"
done

grep -h "$IGNORE_TOKEN" "$TEST_DIR"/* |
    sed 's/^Experiment,//' \
    > ${EXP_DIR}/res_ga_sr.csv

python3 ${EXP_DIR}/notify.py

printf "\nDone.\n"
