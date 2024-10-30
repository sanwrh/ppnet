import argparse
import os
import numpy as np
import pandas as pd
import torch
import pickle  # Import pickle for saving/loading
from Bio import SeqIO
from Bio.PDB import Polypeptide
from transformers import BertTokenizer, BertModel, pipeline
from pepnn_seq.models import FullModel
from pepnn_seq.models import score
from scipy.stats import ttest_ind
import shutil
def determine_final_result(num_interacting_residues, thresholds=range(50), results=range(50)):
    for threshold, result in zip(thresholds, results):
        if num_interacting_residues == threshold:
            return result
    return -1  # 返回 -1 表示超出范围

def load_uniprot_mapping(file_path):
    # 假设该函数实现了从 Excel 文件加载 UniProt 映射的功能
    return pd.read_excel(file_path).set_index('PDBID')['UniProtID'].to_dict()

def analyze_interactions(results_folder, target_folder, uniprot_mapping_file, significance_threshold, weighted_value_threshold,interaction_threshold,analysis_file):
    num = 1

    uniprot_mapping = load_uniprot_mapping(uniprot_mapping_file)
    results = []

    for peptide_folder in os.listdir(results_folder):
        peptide_path = os.path.join(results_folder, peptide_folder)
        if os.path.isdir(peptide_path):
            output_file = os.path.join('example/out/peptides', f"{peptide_folder}_interaction.csv")
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            print(num)
            num += 1

            for seq_folder in os.listdir(peptide_path):
                seq_path = os.path.join(peptide_path, seq_folder)
                binding_site_file = os.path.join(seq_path, 'binding_site_prediction.csv')
                target_file = os.path.join(target_folder, f"{seq_folder.replace('_seq', '')}_site.csv")
                protein_file = os.path.join(target_folder, f"{seq_folder.replace('_seq', '')}_protein.csv")

                if os.path.exists(binding_site_file) and os.path.exists(target_file) and os.path.exists(protein_file):
                    interaction_data = pd.read_csv(binding_site_file)
                    active_sites_data = pd.read_csv(target_file)
                    protein_data = pd.read_csv(protein_file)

                    updated_interaction_data = interaction_data.copy()
                    if len(protein_data) >= len(interaction_data):
                        updated_interaction_data['Position'] = protein_data['Position'].iloc[:len(interaction_data)].values
                    else:
                        raise ValueError("protein_data 行数少于 interaction_data，无法进行匹配")

                    merged_data = pd.merge(active_sites_data, updated_interaction_data, on=['Position', 'Amino acid'], how='inner')

                    if not merged_data.empty:
                        output_folder = os.path.join('example/out/peptides', peptide_folder)
                        os.makedirs(output_folder, exist_ok=True)
                        site_file = os.path.join(output_folder, f"{seq_folder.replace('_seq', '')}_site.csv")
                        interaction_file = os.path.join(output_folder, f"{seq_folder.replace('_seq', '')}_interaction.csv")
                        active_sites_file = os.path.join(output_folder, f"{seq_folder.replace('_seq', '')}_active.csv")

                        updated_interaction_data.to_csv(interaction_file, index=False)
                        active_sites_data.to_csv(active_sites_file, index=False)
                        merged_data.to_csv(site_file, index=False)

                        total_active_sites = merged_data.shape[0]
                        overall_probabilities = updated_interaction_data['Probabilities']
                        active_probabilities = merged_data['Probabilities']

                        t_stat, p_value = ttest_ind(overall_probabilities, active_probabilities, equal_var=False)

                        significance = 1 if p_value < significance_threshold else 0
                        overall_interaction_mean = updated_interaction_data['Probabilities'].mean()
                        active_interaction_mean = merged_data['Probabilities'].mean()
                        weighted_value = np.log2(active_interaction_mean / overall_interaction_mean)

                        dynamic_threshold = active_interaction_mean

                        interacting_residues = merged_data[merged_data['Probabilities'] > dynamic_threshold]
                        num_interacting_residues = interacting_residues.shape[0]
                        interaction_ratio = num_interacting_residues / total_active_sites if total_active_sites > 0 else 0
                        average_probability = interacting_residues['Probabilities'].apply(lambda x: np.log(1 + x)).mean() if num_interacting_residues > 0 else 0

                        final_result = 1 if weighted_value >= weighted_value_threshold and active_interaction_mean >= interaction_threshold and significance == 1 else 0

                        uniprot_id = uniprot_mapping.get(seq_folder[:4], 'Not Found')

                        results.append([
                            seq_folder.replace('_seq', ''),
                            interaction_ratio,
                            overall_interaction_mean,
                            active_interaction_mean,
                            weighted_value,
                            average_probability,
                            significance,
                            final_result,
                            uniprot_id
                        ])

                        merged_data.transpose().to_csv(output_file, mode='a', header=False, index=True)

    if results:
        results_df = pd.DataFrame(results, columns=['Folder Name', 'Interaction Ratio', 'overall_interaction_mean',
                                                    'active_interaction_mean', 'Weighted Value', 'Average Probability',
                                                    'significance', 'Final Result', 'UniProtID'])
        results_df.to_csv(os.path.join('example/out/peptides', f"{peptide_folder}_summary_interaction_results.csv"),
                          index=False, mode='w')
  # 如果文件夹为空
    summarize_uniprot_results('example/out/peptides', analysis_file)

def summarize_uniprot_results(final_site_folder,analysis_file):
    uniprot_stats = {}

    for filename in os.listdir(final_site_folder):
        if filename.endswith('_summary_interaction_results.csv'):
            file_path = os.path.join(final_site_folder, filename)
            df = pd.read_csv(file_path)
            processed_uniprots = set()

            for _, row in df.iterrows():
                uniprot_id = row['UniProtID']
                if uniprot_id not in uniprot_stats:
                    uniprot_stats[uniprot_id] = 0

                if uniprot_id not in processed_uniprots:
                    if df[df['UniProtID'] == uniprot_id]['Final Result'].max() == 1:
                        uniprot_stats[uniprot_id] += 1
                    processed_uniprots.add(uniprot_id)

    result_df = pd.DataFrame(list(uniprot_stats.items()), columns=['UniProtID', 'Final Result Sum'])
    result_df.to_csv(analysis_file, index=False)

    # 确保路径存在
    if os.path.exists('example/out/peptides'):
        # 删除文件夹及其内容
        shutil.rmtree('example/out/peptides')
    print(f"统计完成，结果已保存为{analysis_file}")


def to_var(x):
    if torch.cuda.is_available():
        x = x.cuda()
    return x

# Define a function to load the model
def load_model(params_path):
    model = FullModel(6, 64, 6, 64, 128, 64, dropout=0.2)
    if torch.cuda.is_available():
        model.load_state_dict(torch.load(params_path))
    else:
        model.load_state_dict(torch.load(params_path, map_location='cpu'))
    if torch.cuda.is_available():
        model.cuda()

    model.eval()
    return model

# Define a function to process protein sequences
def process_protein_sequences(protein_folder):
    prot_dict = {}
    protbert_dir = os.path.join(os.path.dirname(__file__), 'pepnn_seq/models/ProtBert-BFD/')

    vocabFilePath = os.path.join(protbert_dir, 'vocab.txt')
    tokenizer = BertTokenizer(vocabFilePath, do_lower_case=False)
    seq_embedding = BertModel.from_pretrained(protbert_dir)

    if torch.cuda.is_available():
        seq_embedding = seq_embedding.cuda()
        seq_embedding = pipeline('feature-extraction', model=seq_embedding, tokenizer=tokenizer, device=0)
    else:
        seq_embedding = pipeline('feature-extraction', model=seq_embedding, tokenizer=tokenizer, device=-1)

    for filename in os.listdir(protein_folder):
        if filename.endswith(".fasta"):
            protein_file = os.path.join(protein_folder, filename)
            records = SeqIO.parse(protein_file, format="fasta")
            prot_sequence = ' '.join(list(records)[0].seq)

            embedding = seq_embedding(prot_sequence)
            embedding = np.array(embedding)

            seq_len = len(prot_sequence.replace(" ", ""))
            start_Idx = 1
            end_Idx = seq_len + 1
            seq_emd = embedding[0][start_Idx:end_Idx]

            prot_seq = to_var(torch.FloatTensor(seq_emd).unsqueeze(0))

            if torch.cuda.is_available():
                prot_seq = prot_seq.cuda()
            # Store processed sequence in the dictionary
            prot_dict[filename] = (prot_seq, prot_sequence)
    return prot_dict

# Define a function to process peptide sequences from a folder
def process_peptide_sequences(peptide_folder):
    peptide_dict = {}
    for filename in os.listdir(peptide_folder):
        if filename.endswith(".fasta"):
            peptide_file = os.path.join(peptide_folder, filename)
            records = SeqIO.parse(peptide_file, format="fasta")
            pep_sequence = str(list(records)[0].seq).replace("X", "")
            pep_sequence = [Polypeptide.d1_to_index[i] for i in pep_sequence]
            peptide_name = os.path.splitext(filename)[0]  # Get the peptide name from filename

            pep_seq = to_var(torch.LongTensor(pep_sequence).unsqueeze(0))

            if torch.cuda.is_available():
                pep_seq = pep_seq.cuda()

            peptide_dict[peptide_name] = pep_seq

    return peptide_dict

# Define a function to calculate scores and output results
def output_results(outputs, prot_sequence, output_directory):
    score_prm = score(outputs)

    with open(os.path.join(output_directory, "prm_score.txt"), 'w') as output_file:
        output_file.writelines("The input protein's score is {0:.2f}".format(score_prm))

    outputs = np.exp(outputs[0])
    amino_acids = []
    probabilities = []
    position = []

    for index, aa in enumerate(prot_sequence.split(" ")):
        probabilities.append(outputs[index, 1])
        amino_acids.append(aa)
        position.append(index + 1)

    output = pd.DataFrame({
        "Position": position,
        "Amino acid": amino_acids,
        "Probabilities": probabilities
    })

    output.to_csv(os.path.join(output_directory, "binding_site_prediction.csv"), index=False)

# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--analysis', action='store_true', help='Run in analysis mode')

    # Common parameters for both modes
    # parser = argparse.ArgumentParser(description="Analyze interactions for peptides.")
    parser.add_argument('-r','--results_folder', type=str, default='example/out', help='Path to the results folder')
    parser.add_argument('-t', '--target_folder', type=str, default='target_min', help='Path to the target folder')
    parser.add_argument('-u','--uniprot_mapping_file', type=str, default='UniProt_PDB_lib.xlsx', help='Path to the UniProt mapping file')
    parser.add_argument('-s','--significance_threshold', type=float, default=0.001, help='Significance threshold for t-test')
    parser.add_argument('-w','--weighted_value_threshold', type=float, default=1.5, help='Threshold for weighted value')
    parser.add_argument('-i','--interaction_value_threshold', type=float, default=0.5, help='Threshold for interaction value')
    parser.add_argument('-of', '--output_file', type=str, default='analysis_output.csv',
                        help='Path to the output file for analysis mode')

    # Parameters for prediction mode
    parser.add_argument('-pf',"-prot_folder", dest="input_protein_folder", required=False, type=str,
                        help="Folder containing fasta files with protein sequences", default="target_min/")
    parser.add_argument('-pp',"-pep_folder", dest="input_peptide_folder", required=False, type=str,
                        help="Folder containing fasta files with peptide sequences", default="example/peptides/")
    parser.add_argument('-o', dest="output_directory", required=False, type=str,
                        help="Base output directory", default="example/out")
    args = parser.parse_args()

    if args.analysis:
        analyze_interactions(
            results_folder=args.results_folder,
            target_folder=args.target_folder,
            uniprot_mapping_file=args.uniprot_mapping_file,
            significance_threshold=args.significance_threshold,
            weighted_value_threshold=args.weighted_value_threshold,
            interaction_threshold = args.interaction_value_threshold,
            analysis_file= args.output_file
        )
    else:

        params_path ='pepnn_seq/params/params.pth'
        # Ensure base output directory exists
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)

        # Define the file path for saving/loading prot_dict
        prot_dict_file = os.path.join(args.output_directory, "prot_dict.pkl")

        # Load prot_dict if it exists
        if os.path.exists(prot_dict_file):
            with open(prot_dict_file, 'rb') as f:
                prot_dict = pickle.load(f)
        else:
            # Process protein sequences
            prot_dict = process_protein_sequences(args.input_protein_folder)

            # Save prot_dict to a file
            with open(prot_dict_file, 'wb') as f:
                pickle.dump(prot_dict, f)

        # Process peptide sequences from the specified folder
        peptide_dict = process_peptide_sequences(args.input_peptide_folder)

        # Load model
        model = load_model(os.path.join(os.path.dirname(__file__), params_path))
        num=1
        # Perform prediction for each peptide sequence against all protein sequences
        for peptide_name, pep_seq in peptide_dict.items():
            print(num)
            num=num+1
            abc =1
            for filename, (prot_seq, prot_sequence) in prot_dict.items():
                print(abc)
                abc =abc+1
                # Extract PDB ID from filename (assuming filename format is {PDBID}.fasta)
                pdb_id = os.path.splitext(filename)[0]  # Get PDB ID from filename (without extension)

                # Create output directory for this peptide and protein PDB ID
                pdb_output_directory = os.path.join(args.output_directory, peptide_name, pdb_id)
                if not os.path.exists(pdb_output_directory):
                    os.makedirs(pdb_output_directory)

                # Model inference
                if torch.cuda.is_available():
                    outputs = model(pep_seq, prot_seq).cpu().detach().numpy()
                else:
                    outputs = model(pep_seq, prot_seq).detach().numpy()

                # Output results
                output_results(outputs, prot_sequence, pdb_output_directory)
