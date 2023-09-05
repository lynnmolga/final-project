import pandas as pd
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns


def load_data(file_path):
    """Load data from a CSV file and return the DataFrame."""
    data = pd.read_csv(file_path)
    return data


def filter_data(data):
    """Filter out genes with NaN values and return the filtered DataFrame."""
    filtered_data = data.dropna(subset=["expression", "mean", "std", "z_score"])
    return filtered_data


def create_gene_color_dict(data):
    """Create a dictionary to map genes to colors based on unique gene names."""
    unique_genes = data["gene"].unique()
    gene_colors = sns.color_palette("tab10", n_colors=len(unique_genes))
    gene_color_dict = {gene: color for gene, color in zip(unique_genes, gene_colors)}
    return gene_color_dict


def perform_tsne(group_data):
    """Perform t-SNE for a given group of genes and return the results."""
    X = group_data[["expression", "mean", "std", "z_score"]].values
    tsne = TSNE(n_components=2, random_state=42, perplexity=5)
    X_tsne = tsne.fit_transform(X)
    return X_tsne


def perform_umap(group_data):
    """Perform UMAP for a given group of genes and return the results."""
    X = group_data[["expression", "mean", "std", "z_score"]].values
    umap_model = umap.UMAP(n_components=2, random_state=42)
    X_umap = umap_model.fit_transform(X)
    return X_umap


def cluster_and_plot(data, gene_color_dict, num_clusters=3):
    """Cluster data using k-means and plot t-SNE and UMAP results with clusters."""
    grouped_data = data.groupby("tumor")  # Group data by tumor type

    for tumor_type, group_data in grouped_data:  # Iterate through each tumor type and its respective data
        tsne_result = perform_tsne(group_data)  # Perform t-SNE for the current tumor type's data
        if tsne_result is None:
            continue

        umap_result = perform_umap(group_data)  # Perform UMAP for the current tumor type's data
        if umap_result is None:
            continue

        kmeans = KMeans(n_clusters=num_clusters, random_state=42)  # Initialize KMeans clustering
        tsne_clusters = kmeans.fit_predict(tsne_result)  # Cluster t-SNE results
        umap_clusters = kmeans.fit_predict(umap_result)  # Cluster UMAP results

        plt.figure(figsize=(10, 5))

        # Plot t-SNE results
        plt.subplot(1, 2, 1)
        for gene in gene_color_dict:
            gene_mask = group_data["gene"] == gene
            color_for_gene = gene_color_dict.get(gene, "black")
            plt.scatter(tsne_result[gene_mask, 0], tsne_result[gene_mask, 1], c=[color_for_gene], label=gene, s=10)
        plt.title(f"t-SNE - Tumor Type: {tumor_type}")
        plt.xlabel("t-SNE Dimension 1")
        plt.ylabel("t-SNE Dimension 2")
        plt.legend()

        # Plot UMAP results
        plt.subplot(1, 2, 2)
        for gene in gene_color_dict:
            gene_mask = group_data["gene"] == gene
            color_for_gene = gene_color_dict.get(gene, "black")
            plt.scatter(umap_result[gene_mask, 0], umap_result[gene_mask, 1], c=[color_for_gene], label=gene, s=10)
        plt.title(f"UMAP - Tumor Type: {tumor_type}")
        plt.xlabel("UMAP Dimension 1")
        plt.ylabel("UMAP Dimension 2")
        plt.legend()

        plt.tight_layout()
        plt.savefig(f'{tumor_type}_umap_tsne_clusters.png')  # Save the plot as an image
        plt.show()  # Display the plot


if __name__ == "__main__":
    file_path = "combined_output.csv"
    data = load_data(file_path)
    filtered_data = filter_data(data)
    gene_color_dict = create_gene_color_dict(filtered_data)
    cluster_and_plot(filtered_data, gene_color_dict)
