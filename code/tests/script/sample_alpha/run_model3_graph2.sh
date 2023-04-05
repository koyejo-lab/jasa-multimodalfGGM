#please change the parameter in the bracket
source activate [`conda environment`]
cd [`path to run_model3`] 


for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
do
    mkdir -p -- "../result_sample_batch/ss_${i}"  
    
    for j in 1 2 3 4 6 9;
    do
        python run_model3_N_alpha.py --p 50   --type 'tridiag2' --s 4 --alpha ${j} --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
        python run_model3_N_alpha.py --p 100  --type 'tridiag2' --s 7 --alpha ${j} --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
        python run_model3_N_alpha.py --p 150  --type 'tridiag2' --s 4 --alpha ${j} --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.01 --lr_a 0.0001 --lr_b 0.01
    done
done
