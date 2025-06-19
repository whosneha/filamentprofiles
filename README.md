# Profile Analysis Tool

This Python tool enables users to extract and analyze **density**, **metallicity**, and other scalar field profiles along cosmological filaments, typically generated from 3D simulation data. It allows for flexible box extraction, annular masking, and multiple modes of profile generation. The tool is designed for use with datasets from simulations like **Horizon-AGN**, **MillenniumTNG**, etc.

##  Features

- Extract density/metallicity profiles along filamentary structures.
- Create 3D sub-boxes from data cubes centered on filaments.
- Apply optional **annular masks** to exclude central regions.
- Multiple profile types supported: mean, median, histogram, gradient.
- Modular structure for easy integration with simulation pipelines.

## File Overview

**`profile_analysis.py`**  
Contains core functions to:
- Extract cube data around filaments.
- Apply annular masks to exclude inner cores.
- Compute scalar profiles (mean, median, histogram, gradient) from input fields.
- Visualize or save profiles depending on user preference.

##  Requirements

You’ll need the following Python libraries:

```bash
pip install numpy matplotlib scipy
```

## Usage

```python
from profile_analysis import analyze_profiles

# Example data inputs
data_cube = np.load("path_to_density_cube.npy")  # 3D array
filament_dict = { ... }  # Node-saddle dictionary from your analysis
box_size_mpc = 2.0
box_radius_mpc = 0.5

# Run analysis
analyze_profiles(
    data_cube,
    filament_dict,
    box_size_mpc,
    box_radius_mpc,
    profile_mode="mean",       # Options: "mean", "median", "hist", "grad"
    use_annular_mask=True,     # Whether to block out center
    plot_results=True,         # Enable matplotlib visualization
    save_path="results/"       # Directory for saving profile files
)
```

## Profile Modes

- `"mean"`: Average scalar values along filament length.
- `"median"`: Robust central tendency profiles.
- `"hist"`: Histogram distributions in each segment.
- `"grad"`: Gradient estimates between segments.

##  Inputs

| Parameter         | Description                                           |
|------------------|-------------------------------------------------------|
| `data_cube`       | 3D NumPy array (e.g., density or metallicity)        |
| `filament_dict`   | Nested dict `{node: {saddle: [segment_coords]}}`     |
| `box_size_mpc`    | Total box side length in Mpc                         |
| `box_radius_mpc`  | Radius used for profile region                       |
| `use_annular_mask`| Excludes central radius (optional)                   |
| `plot_results`    | Plots profiles with matplotlib                       |
| `save_path`       | Directory for saving results (optional)              |

## Output

- `.npy` or `.csv` files of computed profiles
- Optional visual plots (PNG/PDF)
- Printed summary statistics in terminal

## Example Applications

- Compare filament density profiles across MHD vs non-MHD simulations.
- Study metallicity evolution in large-scale structure.
- Benchmark simulation resolution using profile gradients.

## Contributing

Interested in extending functionality (e.g., for temperature, velocity fields, or different masking methods)? Feel free to fork and submit a pull request!

## License

MIT License 
