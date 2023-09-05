import pandas as pd
import matplotlib.pyplot as plt


def load_and_filter_data(file_path):
    """Load and filter data from a CSV file."""
    df = pd.read_csv(file_path, encoding='windows-1252')
    return df


def create_gene_boxplots(df):
    """Create box plots of gene expression for different tumor types."""
    # Separate the DataFrame into three DataFrames based on tumor type
    nn_df = df[df['tumor'] == 'NN']
    bn_df = df[df['tumor'] == 'BN']
    mn_df = df[df['tumor'] == 'MN']

    # Get unique genes in the DataFrame
    genes = df['gene'].unique()

    # Iterate over each gene
    for gene in genes:
        # Create a new figure for each gene
        plt.figure(figsize=(15, 10))

        # Filter the DataFrames for the current gene
        nn_gene_df = nn_df[nn_df['gene'] == gene]
        bn_gene_df = bn_df[bn_df['gene'] == gene]
        mn_gene_df = mn_df[mn_df['gene'] == gene]

        # Create the box plot
        data = [nn_gene_df['expression'], bn_gene_df['expression'], mn_gene_df['expression']]
        labels = ['NN Tumor', 'BN Tumor', 'MN Tumor']
        plt.boxplot(data, labels=labels)

        # Set labels and title
        plt.xlabel('Tumor Type')
        plt.ylabel('Expression')
        plt.title(f'Box Plot of Gene Expression: {gene}')

        # Rotate x-axis labels for better readability (optional)
        plt.xticks(rotation=45)

        # Save the plot as an image file
        plt.savefig(f'{gene}_boxplot.png')

        # Display the plot
        plt.show()


if __name__ == "__main__":
    file_path = 'combined_output.csv'
    df = load_and_filter_data(file_path)
    create_gene_boxplots(df)
