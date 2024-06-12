"""
Function to compare alphafold predictions with experimental data.

Functions:
load_and_process_data
plot_violin_plots
perform_logistic_regression
plot_roc_curves
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
import numpy as np

def load_and_process_data(alphafold_path, experimental_path):
    """
    Load and process data from AlphaFold and experimental datasets.

    Parameters:
        - alphafold_path (str): Path to the AlphaFold data file.
        - experimental_path (str): Path to the experimental data file.

    Returns:
        - pd.DataFrame: Merged dataset with standardized protein pair names and evidence column
          which records whether experimental data is available for the protein pair.
    """
    def load_data(file_path):
        """Helper function to load data from a file, supporting CSV and Excel formats."""
        if file_path.endswith('.xlsx'):
            return pd.read_excel(file_path)
        elif file_path.endswith('.csv'):
            return pd.read_csv(file_path)
        else:
            raise ValueError("Unsupported file format. Please use either a '.csv' or '.xlsx' file.")
    # Load data from files
    alphafold_data = load_data(alphafold_path)
    experimental_data = load_data(experimental_path)

    # Process AlphaFold data
    alphafold_data['Protein_Pair'] = alphafold_data.apply(lambda row: row['Protein1'] + '-' + row['Protein2'], axis=1)
    grouped = alphafold_data.groupby('Protein_Pair').agg({
        'ipTM': 'max',
        'min_PAE': 'min',
        'pDockQ': 'max',
        'ppv': 'max'
    }).reset_index()

    # Process experimental data
    experimental_data = experimental_data.drop_duplicates()
    experimental_data['Protein_Pair'] = experimental_data.apply(lambda row: row['Official Symbol Interactor A'] + '-' + row['Official Symbol Interactor B'], axis=1)

    # Standardize protein pair names
    for dataset in [grouped, experimental_data]:
        dataset['Sorted_Protein_Pair'] = dataset['Protein_Pair'].apply(lambda x: '-'.join(sorted(x.split('-'))))

    # Merge datasets on standardized protein pair names
    merged_data = grouped.merge(experimental_data[['Sorted_Protein_Pair', 'Experimental System']],
                                on='Sorted_Protein_Pair', how='left')
    merged_data['Evidence'] = merged_data['Experimental System'].notna().replace({True: 'Yes', False: 'No'})
    
    return merged_data.drop(columns=['Sorted_Protein_Pair']).drop_duplicates()

def plot_violin_plots(data):
    # Plotting violin plots for various scores
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for i, metric in enumerate(['ipTM', 'pDockQ', 'min_PAE']):
        sns.stripplot(ax=axes[i], x='Evidence', y=metric, data=data, jitter=0.1, color='black', alpha=0.3)
        sns.violinplot(ax=axes[i], x='Evidence', y=metric, data=data, inner=None, scale='width')
        if metric in ['ipTM', 'pDockQ']:
            axes[i].set_ylim(0, 1)
        else:
            axes[i].set_ylim(0, 30)
        axes[i].set_title(f'{metric} Scores by Evidence')
        axes[i].set_xlabel('Evidence')
        axes[i].set_ylabel(f'{metric} Score')
    plt.tight_layout()
    plt.show()

def perform_logistic_regression(data):
    features = ['ipTM', 'pDockQ', 'min_PAE']
    models = {}
    auc_scores = {}
    test_sets = {}
    
    for feature in features:
        # Drop rows with NaN values in the current feature and the 'Evidence' column
        subset_data = data[[feature, 'Evidence']].dropna()

        # Extract X and y for the current feature
        X = subset_data[[feature]]
        y = (subset_data['Evidence'] == 'Yes').astype(int)

        X_train, X_test, y_train, y_test = train_test_split(X[[feature]], y, test_size=0.2, random_state=42)
        model = LogisticRegression().fit(X_train, y_train)
        auc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
        models[feature] = model
        auc_scores[feature] = auc
        test_sets[feature] = (y_test, X_test)
        
    return models, auc_scores, test_sets

def plot_roc_curves(models, auc_scores, test_sets):
    plt.figure(figsize=(12, 6))
    for metric, model in models.items():
        y_test, X_test = test_sets[metric]
        fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
        plt.plot(fpr, tpr, label=f'{metric} (AUC = {auc_scores[metric]:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Different Metrics')
    plt.legend(loc='lower right')
    plt.show()

def plot_roc_curves(models, auc_scores, test_sets):
    plt.figure(figsize=(12, 6))
    for metric, model in models.items():
        y_test, X_test = test_sets[metric]
        probabilities = model.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, probabilities)

        plt.plot(fpr, tpr, label=f'{metric} (AUC = {auc_scores[metric]:.2f})')
        
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Different Metrics')
    plt.legend(loc='lower right')
    plt.show()

#file_path = 'data/Top_predictions.xlsx'
file_path = "data/alphafold_predictions_results.csv"
csv_path = 'data/Merged_remove_duplicate_all.csv'
data = load_and_process_data(file_path, csv_path)
plot_violin_plots(data)
models, auc_scores, test_sets = perform_logistic_regression(data)
print(f"AUC Scores: {auc_scores}")
plot_roc_curves(models, auc_scores, test_sets)
