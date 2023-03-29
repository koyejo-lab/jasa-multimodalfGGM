### Synthesize data

To generate data with isotropic noise (Noise Model 1), run
```
    Rscript datageneration_N1.R cov_name=[graph name] path=[path to store data]
```

Generate the batch of data by running the script file:
```
    chmod a+x run_dgp_N1.sh
    ./run_dgp_N1.sh
```


To generate data with non-isotropic noise (Noise Model 2), run
```
    Rscript datageneration_N2.R cov_name=[graph name] path=[path to store data]
```

Generate the batch of data by running the script file:
```
    chmod a+x run_dgp_N2.sh
    ./run_dgp_N2.sh
```


To generate data with different sample size, run
```
    Rscript datageneration_sample.R cov_name=[graph name] path=[path to store data]
```

Generate the batch of data by running the script file:
```
    chmod a+x run_dgp_sample.sh
    ./run_dgp_sample.sh
```



