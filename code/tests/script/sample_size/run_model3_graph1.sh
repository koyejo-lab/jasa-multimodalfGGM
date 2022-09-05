#please change the parameter in the bracket
source activate [`conda environment`]
cd [`path to run_model3`] 


for i in 1 2 3 4 5 6 7 8 9 10 ;
do
    mkdir -p -- "../result_sample_batch/ss_${i}"  
    python run_model3_N.py           --p 50  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    python run_model3_N_trueparam.py --p 50  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    
    python run_model3_N.py           --p 100  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    python run_model3_N_trueparam.py --p 100  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01

    python run_model3_N.py           --p 150  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
    python run_model3_N_trueparam.py --p 150  --type 'tridiag1' --s 4 --alpha 3 --filepath '../data_sample_batch/ss_${i}'  --id ${i}  --savepath '../result_sample_batch/ss_${i}' --lr_initb 0.1 --lr_a 0.0001 --lr_b 0.01
        
done
