import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import numpy as np

# ────────────────────────────────────────────────
#  Configuration
# ────────────────────────────────────────────────
csv_path = "100_results.csv"           # change if needed
boxplot_output_path = "./output/boxplots_all_params.png"
pca_output_path = "./output/pca_plot.png"
logreg_coeffs_output_path = "./output/logreg_coefficients.png"
logreg_plots_output_path = "./output/logreg_individual_plots.png"

# ────────────────────────────────────────────────
#  Load & prepare data
# ────────────────────────────────────────────────
df = pd.read_csv(csv_path)

# Identify all parameter columns (exclude non-parameters)
exclude_cols = [
    'simulation_index', 'final_time', 'total_cells', 'cell_count_empty',
    'cell_count_susceptible', 'cell_count_infected', 'cell_count_resistant',
    'cell_count_senescent', 'infected_fraction', 'cleared'
]
params_of_interest = [col for col in df.columns if col not in exclude_cols]

# Add 'Outcome' column for visualization
df['Outcome'] = df['cleared'].map({True: 'Cleared', False: 'Persistent'})
palette = {'Cleared': '#2ca02c', 'Persistent': '#d62728'}

# ────────────────────────────────────────────────
#  Boxplots for all parameters
# ────────────────────────────────────────────────
n_params = len(params_of_interest)
n_cols = 5  # Adjust according to taste
n_rows = (n_params + n_cols - 1) // n_cols

fig, axes = plt.subplots(
    nrows=n_rows,
    ncols=n_cols,
    figsize=(n_cols * 3.8, n_rows * 3.2),
    sharex=False,
    sharey=False
)

# Flatten axes for easier iteration
axes_flat = axes.flatten()

for i, param in enumerate(params_of_interest):
    ax = axes_flat[i]
    
    sns.boxplot(
        data=df,
        x='Outcome',
        y=param,
        ax=ax,
        palette=palette,
        width=0.45,
        showfliers=True,
        linewidth=1.4
    )
    
    ax.set_title(param, fontsize=11, pad=8)
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.grid(True, axis='y', linestyle='--', alpha=0.4)

# Hide empty subplots if any
for j in range(i + 1, len(axes_flat)):
    axes_flat[j].set_visible(False)

fig.suptitle("Distribution of All Model Parameters across Monte Carlo Runs\n(boxplots grouped by infection outcome)", 
             fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(boxplot_output_path, dpi=300, bbox_inches='tight')
print(f"Boxplots saved to: {boxplot_output_path}")

# ────────────────────────────────────────────────
#  PCA Plot with low alpha for density visualization
# ────────────────────────────────────────────────
# Standardize parameters
X = df[params_of_interest]
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA to 2 components
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Create PCA dataframe
pca_df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
pca_df['Outcome'] = df['Outcome']

# Plot with low alpha (0.1) to show density, especially for False (Persistent) cases
fig_pca, ax_pca = plt.subplots(figsize=(10, 8))
sns.scatterplot(
    data=pca_df,
    x='PC1',
    y='PC2',
    hue='Outcome',
    palette=palette,
    alpha=0.1,  # Low alpha for visibility of density
    s=50,
    ax=ax_pca
)
ax_pca.set_title("PCA of Parameters Colored by Outcome\n(Low alpha to show density, especially Persistent cases)")
ax_pca.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(pca_output_path, dpi=300, bbox_inches='tight')
print(f"PCA plot saved to: {pca_output_path}")

# ────────────────────────────────────────────────
#  Logistic Regression: Coefficients and Individual Plots
# ────────────────────────────────────────────────
coefficients = {}
logreg_models = {}

for param in params_of_interest:
    # Univariate logistic regression: cleared ~ param
    X_param = df[[param]]
    y = df['cleared'].astype(int)  # Binary: True=1, False=0
    
    model = LogisticRegression()
    model.fit(X_param, y)
    
    # Store coefficient (for the parameter)
    coefficients[param] = model.coef_[0][0]
    
    # Store model for plotting
    logreg_models[param] = model

# Bar plot of coefficients
fig_coeffs, ax_coeffs = plt.subplots(figsize=(12, 8))
coeff_df = pd.DataFrame.from_dict(coefficients, orient='index', columns=['Coefficient'])
coeff_df.sort_values('Coefficient', inplace=True)
coeff_df.plot(kind='barh', ax=ax_coeffs, legend=False)
ax_coeffs.set_title("Logistic Regression Coefficients for Each Parameter\n(Impact on Probability of Cleared=True)")
ax_coeffs.grid(True, axis='x', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(logreg_coeffs_output_path, dpi=300, bbox_inches='tight')
print(f"Logistic regression coefficients saved to: {logreg_coeffs_output_path}")

# Individual logistic plots in a panel
n_rows_log = (n_params + n_cols - 1) // n_cols
fig_log, axes_log = plt.subplots(
    nrows=n_rows_log,
    ncols=n_cols,
    figsize=(n_cols * 4, n_rows_log * 3.5),
    sharex=False,
    sharey=False
)
axes_log_flat = axes_log.flatten()

for i, param in enumerate(params_of_interest):
    ax = axes_log_flat[i]
    
    # Data points
    sns.scatterplot(
        data=df,
        x=param,
        y='cleared'.astype(int),
        ax=ax,
        color='blue',
        alpha=0.3,
        label='Data'
    )
    
    # Logistic curve
    model = logreg_models[param]
    x_vals = np.linspace(df[param].min(), df[param].max(), 300)
    x_vals_2d = x_vals.reshape(-1, 1)
    y_probs = model.predict_proba(x_vals_2d)[:, 1]
    
    ax.plot(x_vals, y_probs, color='red', linewidth=2, label='Logistic Fit')
    
    ax.set_title(f"Prob(Cleared) vs {param}", fontsize=11)
    ax.set_ylabel("Probability of Cleared")
    ax.set_xlabel(param)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.4)

# Hide empty subplots if any
for j in range(i + 1, len(axes_log_flat)):
    axes_log_flat[j].set_visible(False)

fig_log.suptitle("Univariate Logistic Regression Plots for Each Parameter", fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(logreg_plots_output_path, dpi=300, bbox_inches='tight')
print(f"Individual logistic regression plots saved to: {logreg_plots_output_path}")

# Optional: Show plots (comment out if non-interactive)
# plt.show()
