# Code for "Abrupt AMOC threshold behavior in a comprehensive climate model warming run"

Analysis code for the figures in Schulz, Basinski-Ferris, Eisenman and Armour (manuscript submitted to *Nature
Geoscience*). The notebooks start from intermediate data products.

## Running the notebooks

1. Create the environment:

   ```
   conda env create -f environment.yml
   conda activate amoc-threshold
   python -m ipykernel install --user --name amoc-threshold
   ```

2. Put the intermediate data files listed in [data/MANIFEST.csv](data/MANIFEST.csv) into one folder. Data are available
   to reviewers upon request.

3. Set the data folder and start Jupyter from the `notebooks` folder:

   ```
   export AMOC_DATA_DIR=/path/to/data
   export AMOC_FIGURE_DIR=/path/to/figures   # optional; the default is output/ in this folder
   jupyter lab
   ```

Each notebook runs on its own and takes under a minute on the system used here. Figures S6 and S7 read a 3.7 GB file.
Figure S5 and the maps in Figure 2 download the POP grid or map data on first use, so they need internet access.

## Figures

| Figure | Notebook | Input data |
|---|---|---|
| 1 | `fig01_amoc_strength.ipynb` | `amoc_ann_*` |
| 2 | `fig02_climate_impacts.ipynb` | `tas_ann_*`, `pr_mon_*` |
| 3 | `fig03_density_overlap.ipynb` | `FMAASOmean_sigma1_cells_*` |
| 4 | `fig04_density_breakdown_8x.ipynb` | `FMAASOmean_sigma1_cells_8x*`, `FMAmean_TS_cells_8x*`, `surface_TS_sigma1_cells_ctrl*` |
| 5a | `fig05a_surface_freshwater_flux.ipynb` | `SFWF_terms_ann_*_8x` |
| 5b | `fig05b_salinity_section.ipynb` | `SALT_section_atlGIN_8x` |
| S1a | `figS1a_collapsed_amoc.ipynb` | `AMOC_mean_800-2299_*` |
| S1b | `figS1b_subtropical_overlap.ipynb` | `FMAASOmean_sigma1_cells_8x*` |
| S2 | `figS2_density_breakdown_2x4x.ipynb` | `FMAASOmean_sigma1_cells_{2x,4x}*`, `FMAmean_TS_cells_{2x,4x}*`, `surface_TS_sigma1_cells_ctrl*` |
| S3 | `figS3_density_breakdown_60_65N.ipynb` | `ANNmean_TS_cells_8x*`, `surface_TS_sigma1_cells_ctrl*` |
| S4 | `figS4_freshwater_flux_terms.ipynb` | `SFWF_terms_ann_*_{8x,control}` |
| S5 | `figS5_sss_change.ipynb` | `SSS_mean_change_800_2300_8x` |
| S6 | `figS6_isopycnal_streamfunction.ipynb` | `PSI_SIGMA1_POP_MOC_8x` |
| S7 | `figS7_amoc_density_vs_depth.ipynb` | `PSI_SIGMA1_POP_MOC_8x`, `amoc_ann_*8x*` |

Shared settings and calculations are in `notebooks/common.py`: the filter, the area-weighted percentiles, the density
calculation and the regions.

## License

MIT; see [LICENSE](LICENSE).

## Use of AI tools

This code was reformatted for data availability by Claude (Anthropic), which was also used when developing some of the
analysis code, as stated in the manuscript. The authors are fully responsible for the content.