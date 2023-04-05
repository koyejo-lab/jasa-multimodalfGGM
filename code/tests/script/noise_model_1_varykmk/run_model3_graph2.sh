source activate [`conda environment`]
cd [`path to run_model3`] 



for i in 1 2 3 4 5 6 7 8 9 10;
do
    for j in 5 7 9 11;
    do
        for ((k=3;k<=j;k+=2));
        do
            python run_model3_kmk.py --init 'est' --p 50  --N 100  --k ${k} --km ${j} --type 'tridiag2' --filepath "../data_kmk_batch/ss_${i}"  --savepath "../results_kmk_batch/ss_${i}" --noise '' --thre 0.01 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01  
            python run_model3_kmk.py --init 'est' --p 100 --N 100  --k ${k} --km ${j} --type 'tridiag2' --filepath "../data_kmk_batch/ss_${i}"  --savepath "../results_kmk_batch/ss_${i}" --noise '' --thre 0.01 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01 
            python run_model3_kmk.py --init 'est' --p 150 --N 100  --k ${k} --km ${j} --type 'tridiag2' --filepath "../data_kmk_batch/ss_${i}"  --savepath "../results_kmk_batch/ss_${i}" --noise '' --thre 0.01 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01 
        done
    done

done