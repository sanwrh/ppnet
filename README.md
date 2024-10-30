# Network Pharmacology Tools Based on PEPNN

The network pharmacology tools used for predicting peptide-target interactions are based on PEPNN. This script utilizes the PEPNN protein-peptide interaction model as a foundation and constructs a network pharmacology tool using drug target data from the scPDB database. This tool can predict interactions with 5,715 PDB data and 898 protein targets.

You need to install PEPNN according to the guidelines provided at [PepNN Installation](https://gitlab.com/oabdin/pepnn). After installation, proceed with predicting protein-peptide interactions.

Once the PEPNN data runs smoothly, copy all files to the main PepNN folder.

## Modes of ppnet.py

- `python ppnet.py` is for prediction mode.
  - `-pf`, `-prot_folder`: Input folder for the protein target library, with two options: `target` contains all 5,715 PDB data, and `target_mini` contains the 898 PDB data with the lowest resolution. The default is `target_min/`.
  - `-pp`, `-pep_folder`: Folder containing FASTA files of peptide sequences. The default is `example/peptides/`.
  - `-o`, `--output_directory`: Folder for PepNN prediction results, defaulting to `example/out`.
  
**Note**: The first time you run the prediction mode, a file named `prot_dict.pkl` will be generated in the folder set by `-o`. This may take some time. When you run it a second time, it will directly read the `prot_dict.pkl` file from the `-o` directory. If you change the input protein target library (e.g., from `target_min/` to `target`), you need to delete the `prot_dict.pkl` file in the `-o` directory.

- `python ppnet.py -a` is for analysis mode.
  - `-r`, `--results_folder`: The folder for result data, the same as the `-o` folder from the previous step. The default is `example/out`.
  - `-t`, `--target_folder`: The protein target library folder from the previous step. The default is `target_min/`.
  - `-u`, `--uniprot_mapping_file`: The mapping file between UniProt ID and PDB ID. The default is `UniProt_PDB_lib.xlsx`.
  - `-of`, `--output_file`: The output statistical results file. The default is `analysis_output.csv`.
  - `-s`, `--significance_threshold`: The default is 0.001.
  - `-w`, `--weighted_value_threshold`: The default is 1.5.
  - `-i`, `--interaction_value_threshold`: The default is 0.5.

## Key Formulas

1. **T-test**  
   \[
   t_{stat}, p_{value} = ttest\_ind(overall\_probabilities, active\_probabilities, equal\_var=False)
   \]

2. **Significance**  
   \[
   significance = 
   \begin{cases}
   1 & \text{if } p_{value} < significance\_threshold \\
   0 & \text{otherwise}
   \end{cases}
   \]

3. **Overall and Active Interaction Mean**  
   \[
   overall\_interaction\_mean = \text{mean}(updated\_interaction\_data['Probabilities'])
   \]  
   \[
   active\_interaction\_mean = \text{mean}(merged\_data['Probabilities'])
   \]

4. **Weighted Value**  
   \[
   weighted\_value = \log_2\left(\frac{active\_interaction\_mean}{overall\_interaction\_mean}\right)
   \]

5. **Dynamic Threshold**  
   \[
   dynamic\_threshold = active\_interaction\_mean
   \]

6. **Interacting Residues**  
   \[
   interacting\_residues = merged\_data[merged\_data['Probabilities'] > dynamic\_threshold]
   \]  
   \[
   num\_interacting\_residues = \text{shape}(interacting\_residues)[0]
   \]

7. **Interaction Ratio**  
   \[
   interaction\_ratio = \frac{num\_interacting\_residues}{total\_active\_sites} \quad \text{if } total\_active\_sites > 0
   \]

8. **Average Probability**  
   \[
   average\_probability = 
   \begin{cases}
   \text{mean}\left(\log(1 + x)\right) & \text{for } x \in interacting\_residues['Probabilities'] \\
   0 & \text{if } num\_interacting\_residues = 0
   \end{cases}
   \]

9. **Final Result**  
   \[
   final\_result = 
   \begin{cases}
   1 & \text{if } weighted\_value \geq weighted\_value\_threshold \text{ and } active\_interaction\_mean \geq interaction\_threshold \text{ and } significance = 1 \\
   0 & \text{otherwise}
   \end{cases}
   \]

10. **UniProt ID Mapping**  
    \[
    uniprot\_id = uniprot\_mapping.get(seq\_folder[:4], 'Not Found')
    \]
