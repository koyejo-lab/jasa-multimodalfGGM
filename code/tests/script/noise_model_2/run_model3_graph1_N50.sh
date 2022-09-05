#please change the parameter in the bracket
source activate [`conda environment`]
cd [`path to run_model3`] 


for i in   1 2 3 4 5 6 7 8 9 10 ;
do
    mkdir -p -- "../result_N2_batch/ss_${i}"  
    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    


    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    

    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    


    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    
    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 100 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 150 --N 50  --type 'tridiag1' --filepath "../data_N2_batch/ss_${i}"  --savepath "../result_N2_batch/ss_${i}" --noise '_noise3' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001


done

