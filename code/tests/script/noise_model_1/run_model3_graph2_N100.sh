


for i in   1 2 3 4 5 6 7 8 9 10 ;
do
    mkdir -p -- "../result_N1_batch/ss_${i}"  
    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.001 --lr_b 0.01
    


    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.001
    

    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.01 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.001 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.001 --lr_b 0.01
    


    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.001 --lr_a 0.0001 --lr_b 0.01
    
    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    

    python run_model3.py --init 'est' --p 50  --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 100 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001
    python run_model3.py --init 'est' --p 150 --N 100  --type 'tridiag2' --filepath "../data_N1_batch/ss_${i}"  --savepath "../result_N1_batch/ss_${i}" --noise '' --thre 0.001 --lr_initb 0.0001 --lr_a 0.0001 --lr_b 0.001


done

