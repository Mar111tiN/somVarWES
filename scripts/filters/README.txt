filters can be in either R or python
they receive the input table as pd.dataframe and an output file name
- have to be callable as script(input,output)
- have to be executable
- have to be callable as 
- must output (among others) the output file exactly as provided
- parameters of script go into the respective list entry in config['filter'] under params
- params can be updated by directly calling the filter {script} -param1 value input output


! filter with name NAME has to output (beside others) a file with the name output.csv
in order to fulfill output criteria of snakemake rule 


for R filters, 
