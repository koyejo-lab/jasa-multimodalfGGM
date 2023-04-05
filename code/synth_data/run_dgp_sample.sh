for i in  1 2 3 4 5 6 7 8 9 10;
do
    mkdir -p -- "../../data/data_sample_batch/ss_${i}"  
    Rscript datageneration_sample.R cov_name="tridiag1" path="../../data/data_sample_batch/ss_${i}"
    Rscript datageneration_sample.R cov_name="tridiag2" path="../../data/data_sample_batch/ss_${i}"
    Rscript datageneration_sample.R cov_name="tridiag3" path="../../data/data_sample_batch/ss_${i}"

done
