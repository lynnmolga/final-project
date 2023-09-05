import copy
import numpy as np
import pandas as pd

matching_genes_df = pd.DataFrame()

NAMES = {
    ("ny-eso1", "CTAG1B"): "NM_001327",
    "alk": "NM_004304",
    "kras": "NM_033360",
    "egfr": "NM_001346941",
    ("pd-l1", "cd274"): "NM_001267706",
    "pten": "AF067844",
    "tp53": "MH011443",
    "pik3ca": "MT991405",
    "rb1": "L41914",
    "dlec1": ("AB026898", "AB010443"),
    "hopx": "NM_001145460",
    "braf": "NM_001378474",
    "irf1": "NM_002198",
    ("cyfra21", "krt19"): "NM_002276",
    "mir550a3": "NR_039600",
    "mir629": None,
    "ros1": (
        "AH002964",
        "M13368",
        "M13591",
        "M13592",
        "M13593",
        "M13594",
        "M13595",
        "M13596",
        "M13597",
        "M13598",
        "M13599"
    ),
}


def search_name(name, df):
    """
    Search for a gene name or accession in a DataFrame and return the corresponding Probe_Id.

    Args:
        name: Gene name or accession to search for.
        df: The DataFrame containing gene information.

    Returns:
        str, list, or None: If found, returns the Probe_Id as a string or list.
                           If not found, returns None.
    """

    matching_genes = pd.DataFrame()  # Initialize an empty DataFrame to store matching genes
    accession = []  # Initialize an empty list to store accession results

    # Search for gene name(s) in the 'ILMN_Gene' column
    if isinstance(name, tuple):
        for n in name:
            temp = df[df['ILMN_Gene'] == n]  # Search for exact gene name match
            if not temp.empty:
                matching_genes = pd.concat([matching_genes, temp])
            else:
                temp = df[df['ILMN_Gene'] == n.upper()]  # Search for uppercase gene name match
                if not temp.empty:
                    matching_genes = pd.concat([matching_genes, temp])
    else:
        temp = df[df['ILMN_Gene'] == name.upper()]  # Search for uppercase gene name match
        if not temp.empty:
            matching_genes = pd.concat([matching_genes, temp])

    # If matching genes were found, return their 'Probe_Id'(s)
    if len(matching_genes) > 0:
        return matching_genes['Probe_Id'].values
    else:
        # If not found, check for accession match
        if isinstance(name, tuple):
            for n in name:
                accession = df[df['Accession'] == n]  # Search for accession match
            if len(accession) > 0:
                return accession['Probe_Id'].values[0]  # Return the first 'Probe_Id' for the matched accession
        else:
            accession = df[df['Accession'] == name]  # Search for accession match
            if len(accession) > 0:
                return accession['Probe_Id'].values[0]  # Return the 'Probe_Id' for the matched accession

    # If no matches were found, return None
    return None


def search_probe(probe_id):
    """
    Search for rows containing a specific probe_id in a CSV file and store them in a global DataFrame.

    Args:
        probe_id: The probe_id to search for in the CSV file.

    """
    global matching_genes_df
    file_path = "GSE135304_raw.csv"  # Update with the actual file path to the CSV file
    df = pd.read_csv(file_path)  # Read the CSV file into a DataFrame

    # Extract rows with ILMN_0000 number in the first column
    matching_rows = df[df.iloc[:, 0].str.contains(probe_id)]

    if len(matching_genes_df) > 0:
        # Append the matching rows to the global matching_genes_df
        matching_genes_df = pd.concat([matching_genes_df, matching_rows])
    else:
        # Store the matching rows in a new DataFrame if it's empty
        matching_genes_df = pd.DataFrame(matching_rows.values, columns=df.columns)


class Info:
    def __init__(self, tumor, GSM_num, expression, probe_num, flag=False):
        self.type_tumor = tumor
        self.GSM_id = GSM_num
        self.gene_expression = expression
        self.is_pval = flag
        self.probe = probe_num

    def calculate_std(self, objects):
        """
        Calculate standard deviation and z-scores for gene expression values.

        Args:
            objects (list of Info): A list of Info objects with gene expression data.

        Returns:
            list of dict: A list of dictionaries containing probe, GSM, tumor type, expression, mean, std, and z-score.
        """
        std_data = []  # Initialize a list to store calculated standard deviation data
        for obj in objects:
            expressions = [o.gene_expression for o in objects if o.probe == obj.probe]
            # Check if there are less than 2 values for the same probe
            if len(expressions) < 2:
                std_data.append({'probe': obj.probe, 'gsm': obj.GSM_id, 'expression': obj.gene_expression, 'std': None})
            else:
                mean = sum(expressions) / len(expressions)  # Calculate the mean
                variance = sum((x - mean) ** 2 for x in expressions) / (len(expressions) - 1)  # Calculate variance
                std = variance ** 0.5  # Calculate standard deviation
                z_scores = (obj.gene_expression - mean) / std  # Calculate z-scores
                std_data.append({'probe': obj.probe, 'gsm': obj.GSM_id, 'tumor': obj.type_tumor,
                                 'expression': obj.gene_expression, 'mean': mean, 'std': std, 'z_score': z_scores})
        return std_data  # Return the list of calculated data


def create_gene_dataframe(df):
    """
    Create a DataFrame containing gene information based on gene names and NM numbers.

    Args:
        df: The DataFrame containing gene information.

    Returns:
        DataFrame: A DataFrame with gene names as columns and corresponding probe_ids as values.
    """
    gene_data = {}  # Initialize a dictionary to store gene data

    for gene_names, nm_numbers in NAMES.items():
        if isinstance(gene_names, str):
            gene_names = (gene_names,)  # Convert single gene name to tuple for consistency
        if nm_numbers is not None:
            probe_ids = search_name(gene_names, df)  # Search for probe_ids based on gene names
            if probe_ids is not None:
                for gene_name in gene_names:
                    if gene_name not in gene_data:
                        gene_data[gene_name] = []  # Initialize an empty list for each gene name
                    gene_data[gene_name].extend(probe_ids)  # Extend the list with probe_ids

    gene_df = pd.DataFrame.from_dict(gene_data, orient='index').T  # Create a DataFrame from the gene data dictionary
    return gene_df  # Return the gene DataFrame


def main():
    # Read the CSV file into a dataframe
    file_path = "GPL10558_HumanHT-12_V4_0_R1_15002873_B.csv"  # Update with the actual file path
    df = pd.read_csv(file_path)
    # Example search for a specific name
    for name, nm_numbers in NAMES.items():
        results = search_name(name, df)
        if results is not None:
            for r in results:
                search_probe(r)
            print("Results found for:", name, "-")
            print(results)
        else:
            print("No results found for:", name)
    # Filter rows with p-value lower than 0.05
    for i in range(2, len(matching_genes_df.columns), 2):
        pval_col = matching_genes_df.columns[i]
        prev_col = matching_genes_df.columns[i - 1]

        matching_genes_df.loc[matching_genes_df[pval_col] > 0.05, pval_col] = None
        matching_genes_df.loc[matching_genes_df[pval_col] > 0.05, prev_col] = None
    #print(matching_genes_df)
    info_arr = []
    gsm_only = "GSM"
    flag = False
    file_path2 = "series_matrix.csv"  # Update with the actual file path
    df2 = pd.read_csv(file_path2, encoding='windows-1252')
    # for i in range(len(matching_genes_df.index)):
    for i in range(matching_genes_df.shape[0]):  # Iterate over rows
        for j in range(1, matching_genes_df.shape[1], 2):
            if matching_genes_df.iloc[i, j+1] is not None:
                expression_val = matching_genes_df.iloc[i, j]
                expression = expression_val
                flag = True
                GSM_num = int((j + 1)/2 + 4003136)
                gsm_string = gsm_only + str(GSM_num)
                row_index, col_index = np.where(df2 == gsm_string)

                # Retrieve the value at the modified row and col_index
                tumor = df2.iloc[10, col_index].item().split(": ")[1]
                probe = matching_genes_df.iloc[i, 0]
                # Create an Info object with the extracted information
                info_obj = Info(tumor, gsm_string, expression, probe, flag)
                # Append the object to the matrix
                info_arr.append(info_obj)
            else:
                info_arr.append(None)

    gene_df = create_gene_dataframe(df)
    sample_index = 0
    # Create a 712x18 matrix filled with zeros
    result_matrix = np.zeros((713, len(gene_df.columns)), dtype=int)

    # Iterate through matching_genes_df (which is now data) starting from column 1 and jumping 2 columns at a time
    for sample_num in range(1, len(matching_genes_df.columns), 2):
        sample_id_ref = matching_genes_df.columns[sample_num]
        sample_detection_pval = matching_genes_df.columns[sample_num + 1]
        sample_index += 1

        for gene_index, (_, probes) in enumerate(gene_df.items()):
            for probe in probes:
                if probe is not None:
                    probe_found = matching_genes_df["ID_REF"].str.contains(probe).any()

                    if probe_found:
                        # Save the row number where the probe was found
                        probe_row_number = matching_genes_df[matching_genes_df["ID_REF"].str.contains(probe)].index[0]
                    try:
                        # Check if the detection_pval is not None (2 columns down from the sample_num column)
                        if not pd.isna(matching_genes_df[sample_detection_pval][probe_row_number]):
                            # Set the corresponding element in the result_matrix to 1
                            result_matrix[sample_index, gene_index] = 1
                    except Exception as e:
                        # Handle any exceptions here (e.g., index out of bounds, column not found)
                        print(f"Error: {e}")
    # Print the result matrix
    print(result_matrix)

    # Create lists to store Info objects for different tumor types
    BN_arr = []
    MN_arr = []
    NN_arr = []

    # Iterate through info_arr and separate Info objects by tumor type
    for cell in info_arr:
        if cell is not None:
            if cell.type_tumor == "BN":
                cell_copy = copy.deepcopy(cell)
                BN_arr.append(cell_copy)
            if cell.type_tumor == "MN":
                cell_copy = copy.deepcopy(cell)
                MN_arr.append(cell_copy)
            if cell.type_tumor == "NN":
                cell_copy = copy.deepcopy(cell)
                NN_arr.append(cell_copy)

    # Create a DataFrame to store the gene expression data per patient
    gene_per_patient = pd.DataFrame(result_matrix)
    gene_per_patient['tumor type'] = None

    # Populate the 'tumor type' column of gene_per_patient DataFrame
    for cell in info_arr:
        if cell is not None:
            gsm_id = int(cell.GSM_id[3:]) - 4003136
            tumor_local = cell.type_tumor
            gene_per_patient.at[gsm_id, 'tumor type'] = tumor_local

    # Save the gene_per_patient DataFrame to a CSV file
    gene_per_patient.to_csv('genes_per_patient1.csv', index=False)

    # Calculate standard deviations and z-scores for each tumor type
    std_BN = Info.calculate_std(BN_arr[0], BN_arr)
    std_MN = Info.calculate_std(MN_arr[0], MN_arr)
    std_NN = Info.calculate_std(NN_arr[0], NN_arr)

    # Combine the standard deviation data for all tumor types
    all_std_data = []
    all_std_data.extend(std_BN)
    all_std_data.extend(std_MN)
    all_std_data.extend(std_NN)

    # Create a DataFrame from the combined std_data
    combined_df = pd.DataFrame(all_std_data)

    # Save the DataFrame to a CSV file
    combined_df.to_csv('combined_output.csv', index=False)


if __name__ == '__main__':
    main()

