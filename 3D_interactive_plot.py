import plotly.graph_objects as go
import pandas as pd
import plotly.io as pio
# Set the default renderer to browser if running as a standalone script
pio.renderers.default = 'browser'

def plot_3d_interactive(data, x_var, y_var, z_var, color_var, title='3D Scatter Plot'):
    """
    Creates an interactive 3D scatter plot with hovering functionality.

    Parameters:
        data (DataFrame): The pandas DataFrame containing the data.
        x_var (str): Column name for x-axis values.
        y_var (str): Column name for y-axis values.
        z_var (str): Column name for z-axis values.
        color_var (str): Column name for coloring the markers.
        title (str): Title of the plot.
    """
    # Viridis color scale inverted to make high values more visible
    custom_color_scale = [
        [0.0, '#00CED1'],  # Greeny blue (Intense teal)
        [0.25, '#ADFF2F'], # Yellowy green (Bright lime green)
        [0.5, '#FFBF00'],  # Deep orangey yellow (Vibrant amber)
        #[0.75, '#FFA500'], # Orange (Bright orange)
        [0.75, '#FF4500'], # Orange (Bright orange)
        [1.0, '#FF0000']   # Red (Vivid red)
        #[1.0, '#800000']   # Red (Vivid red)
    ]

    # Create a 3D scatter plot
    fig = go.Figure(data=[go.Scatter3d(
        x=data[x_var],
        y=data[y_var],
        z=data[z_var],
        mode='markers',
        marker=dict(
            size=4,
            color=data[color_var],  # Use color_var for the color scale
            colorscale=custom_color_scale,  # Color scale can be adjusted
            colorbar=dict(title=color_var),
            opacity=0.5
        ),
        text=[f'{p1},{p2}' for p1, p2 in zip(data['Protein1_Domain'], data['Protein2_Domain'])],
        hoverinfo='text'
    )])

    # Update the layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title=x_var,
            yaxis_title=y_var,
            zaxis_title=z_var
        )
    )

    # Show the plot
    fig.show()

def plot_3d_interactive_by_category(data, x_var, y_var, z_var, color_var, category_var, title='3D Scatter Plot'):
    fig = go.Figure()

    custom_color_scale = [
        [0.0, '#00CED1'],  # Greeny blue (Intense teal)
        [0.25, '#ADFF2F'], # Yellowy green (Bright lime green)
        [0.5, '#FFBF00'],  # Deep orangey yellow (Vibrant amber)
        #[0.75, '#FFA500'], # Orange (Bright orange)
        [0.75, '#FF4500'], # Orange (Bright orange)
        [1.0, '#FF0000']   # Red (Vivid red)
        #[1.0, '#800000']   # Red (Vivid red)
    ]

    # Filter data by category for separate traces
    for category in data[category_var].unique():
        cat_data = data[data[category_var] == category]
        if category == "Y":
            marker_style = dict(size=10, symbol='circle')
            edge_color = 'rgba(0,0,0,2)'  # Black border for category A
        else:
            marker_style = dict(size=10, symbol='diamond')
            edge_color = 'rgba(255,255,255,2)'  # White border for category B
        
        fig.add_trace(go.Scatter3d(
            x=cat_data[x_var],
            y=cat_data[y_var],
            z=cat_data[z_var],
            mode='markers',
            marker=dict(
                color=cat_data[color_var],
                colorscale=custom_color_scale,
                line=dict(color=edge_color, width=2),  # Edge color based on category
                **marker_style
            ),
            text=[f'Protein1 Domain: {p1}, Protein2 Domain: {p2}' for p1, p2 in zip(cat_data['Protein1_Domain'], cat_data['Protein2_Domain'])],
            hoverinfo='text',
            name=f'Category {category}'  # Label trace by category
        ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title=x_var,
            yaxis_title=y_var,
            zaxis_title=z_var
        )
    )

    fig.show()

# Example Usage
# Load data
data = pd.read_csv('data/alphafold_predictions_results.csv')

# Variables to plot
x_var = 'ipTM'
y_var = 'min_PAE'
z_var = 'pDockQ'
color_var = 'Num_Consistent'
category_var = 'Column1'

# Call the function
#filter data to remove min_PAE >15
data_minPAE_15 = data[data['min_PAE'] <= 15]
plot_3d_interactive(data_minPAE_15, x_var, y_var, z_var, color_var, title='Protein Interaction Plot')

data_minPAE_15 = data[data['min_PAE'] <= 15]
plot_3d_interactive_by_category(data_minPAE_15, x_var, y_var, z_var, color_var, category_var, title='Protein Interaction Plot')
