for i in  1 2 3 4 5 6 7 8 9 10;
do
    mkdir -p -- "../data_N1_batch/ss_${i}"
    Rscript datageneration_N1.R cov_name="tridiag1" path="../data_N1_batch/ss_${i}"
    Rscript datageneration_N1.R cov_name="tridiag2" path="../data_N1_batch/ss_${i}"
    Rscript datageneration_N1.R cov_name="tridiag3" path="../data_N1_batch/ss_${i}"
    Rscript datageneration_N1.R cov_name="power"    path="../data_N1_batch/ss_${i}"

done
