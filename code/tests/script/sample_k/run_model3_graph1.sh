source activate [`conda environment`]
cd [`path to run_model3`] 


for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
do

    python run_model3_N_varyk.py --p 50  --N 324  --type 'tridiag1' --s 4 --alpha 3 --filepath "../data_k_batch/ss_${i}"  --id $i  --savepath '../result_k_batch/' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    python run_model3_N_varyk.py --p 100 --N 324  --type 'tridiag1' --s 4 --alpha 3 --filepath "../data_k_batch/ss_${i}"  --id $i  --savepath '../result_k_batch/' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    python run_model3_N_varyk.py --p 150 --N 289  --type 'tridiag1' --s 4 --alpha 3 --filepath "../data_k_batch/ss_${i}"  --id $i  --savepath '../result_k_batch/' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01

done