#  For admixture we have to convert the LD_pruned data for ADMIXTURE compatibility
#Convert chromosome column to integer
for dataset in gardenonly_pruned
do
    awk '{$1="1"; print}' ${dataset}.bim > tmp.bim
    mv tmp.bim ${dataset}.bim
done


# Run admixture with 10 replicates per K, for K=1 to 10, on LD-pruned dataset

mkdir -p admixture_pruned_replicates
cd admixture_pruned_replicates

for K in 1 2 3 4 5 6 
do
    for rep in 1 2 3 4 5 6 7 8 9 10
    do
        echo "Running K=${K}, replicate=${rep}..."
        admixture --cv -s ${rep} ../gardenonly_pruned.bed $K > log_K${K}_rep${rep}.txt
        
        mv gardenonly_pruned.${K}.Q gardenonly_pruned_K${K}_rep${rep}.Q
        mv gardenonly_pruned.${K}.P gardenonly_pruned_K${K}_rep${rep}.P
    done
    done

    #Extract all CV errors into one file
cd admixture_pruned_replicates

for file in log_K*.txt
do
    K=$(echo $file | sed 's/log_K\([0-9]*\)_rep.*/\1/')
    rep=$(echo $file | sed 's/.*_rep\([0-9]*\)\.txt/\1/')
    cv=$(grep "CV error" $file | sed 's/.*: //')
    echo "$K $rep $cv"
done > all_cv_errors.txt

#we compared our k value the pick our best fit data
#Eventually the final visualization were K=4 ("gardenonly_pruned_clean_maf01_K5_rep3.Q"), K=5 ("gardenonly_pruned_clean_maf01_K5_rep3.Q") and K=6 ("gardenonly_pruned_clean_maf01_K6_rep7.Q")