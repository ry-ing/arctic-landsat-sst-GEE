#!/usr/bin/env python
"""
Derives SST retreival coefficients (A, B, C, D) from libRadtran simulations
and validates against EasyCORA in-situ matchups.

Algorithm: SST = A·Tb_scaled + B·(Tb_scaled·TCWV_scaled) + C + D
Tb and TCWV are independently RobustScaler-normalised before the interaction
term is computed; D is derived
via k-fold (5-fold) cross-validation (see https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html)
"""

# import modules
import pandas as pd
import numpy as np
from pathlib import Path
import glob
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import RobustScaler
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
import joblib
import logging
import sys
import warnings
import json

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from scipy import stats
import seaborn as sns

import geopandas as gpd
from shapely.geometry import Point

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ============================================================================
# 1. Configuration
# ============================================================================
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


# --- Inputs ---

REPO_ROOT = Path(__file__).parent.parent # back one level from the script location
SIMULATION_RESULTS_CSV  = REPO_ROOT / 'data/libradtran_sims/libradtran_simulations_df_MERRA.csv'
VALIDATION_DATA_DIR     = REPO_ROOT / 'data/GEE_validation_files'
COASTLINE_SHP           = REPO_ROOT / 'data/coastline_shapefile/ne_10m_land.shp'
COUNTRY_SHP             = REPO_ROOT / 'data/coastline_shapefile/ne_110m_admin_0_countries.shp'


# --- Outputs ---
(REPO_ROOT / 'logs').mkdir(exist_ok=True)
(REPO_ROOT / 'Figures').mkdir(exist_ok=True)

DIAGNOSTIC_PLOT_FILE    = str(REPO_ROOT / 'Figures/diagnostic_histograms.png')
TRAINING_PLOT_FILE      = str(REPO_ROOT / 'Figures/training_performance_scatter.png')
VALIDATION_PLOT_FILE    = str(REPO_ROOT / 'Figures/validation_scatter_plots.png')
VALIDATION_PLOT_GREENLAND = str(REPO_ROOT / 'Figures/validation_scatter_greenland.png')
VALIDATION_MAP_FILE     = str(REPO_ROOT / 'Figures/validation_location_histograms.png')
VALIDATION_BOXPLOT_FILE        = str(REPO_ROOT / 'Figures/validation_residuals_boxplot.png')
MASTER_PARAMS_CSV       = str(REPO_ROOT / 'data/algorithm_params/MASTER_ALGORITHM_PARAMETERS.csv')
VALIDATION_STATS_CSV    = str(REPO_ROOT / 'data/algorithm_params/VALIDATION_STATISTICS.csv')

# --- Algorithm constants ---
SENSORS_TO_PROCESS = ['L4', 'L5', 'L7', 'L8', 'L9']
SST_MIN_KELVIN = 272.0    # exclude likely sea-ice profiles from MERRA-2 (skin T tail below ~272 K)

# --- Contamination filter ---
BT_STDEV_THRESHOLD = 1.0  # K — max within-scene BT stdev

# --- Cross-validation settings for D ---
CV_N_SPLITS  = 5          # number of k-fold splits
CV_N_SIGMA   = 1          # sigma threshold for outlier removal within each calibration fold
CV_RANDOM_STATE = 42

# order of mag. format colourbar tick labels in validation scatter plot
class OOMFormatter(mticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mticker.ScalarFormatter.__init__(self, useOffset=offset, useMathText=mathText)
    def _set_order_of_magnitude(self): self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(str(REPO_ROOT / 'logs/derivation_validation.log')),
            logging.StreamHandler(sys.stdout),
        ],
    )

def get_sensor_from_id(landsat_id: str) -> str:
    if pd.isna(landsat_id): return None
    if landsat_id.startswith("LC09"): return "L9"
    if landsat_id.startswith("LC08"): return "L8"
    if landsat_id.startswith("LE07"): return "L7"
    if landsat_id.startswith("LT05"): return "L5"
    if landsat_id.startswith("LT04"): return "L4"
    return None

def compute_validation_stats(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Compute RMSE, median accuracy, and median precision."""
    residuals = y_pred - y_true
    rmse      = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    accuracy  = float(np.median(np.abs(residuals)))
    precision = float(np.median(np.abs(residuals - np.median(residuals))))
    bias      = float(np.mean(residuals))
    return {'rmse': rmse, 'accuracy': accuracy, 'precision': precision, 'bias': bias}

# ============================================================================
# 1. Derive A, B, C from radiative transfer simulations
# ============================================================================

def derive_coefficients(simulation_df: pd.DataFrame):
    """Fit A, B, C from libRadtran simulations using an interaction-term linear model."""
    logging.info("--- PHASE 1: Deriving A, B, C from radiative transfer simulations ---")

    master_params = {}
    scalers = {}
    training_df = simulation_df.copy()
    tcwv_col = 'tcwv_merra2_kg_m2'
    y = simulation_df['sst_skin_merra2_k']

    for sensor in SENSORS_TO_PROCESS:
        bt_col = f'sim_bt_{sensor}'
        if bt_col not in simulation_df.columns or tcwv_col not in simulation_df.columns:
            logging.warning(f"  Skipping {sensor}: missing columns.")
            continue

        # Feature columns for scaling (INDEPENDANT scalers)
        feature_cols = [bt_col, tcwv_col]
        df_s = pd.concat([y, simulation_df[feature_cols]], axis=1).dropna()
        if df_s.empty:
            continue

        y_s = df_s['sst_skin_merra2_k']
        X_s = df_s[feature_cols]

        # Fit RobustScaler on the base features
        scaler = RobustScaler()
        X_scaled = scaler.fit_transform(X_s.values)
        scalers[sensor] = scaler

        # Create interaction term: Tb_scaled * TCWV_scaled
        tb_scaled = X_scaled[:, 0]
        tcwv_scaled = X_scaled[:, 1]
        interaction = tb_scaled * tcwv_scaled

        # Fit model with 3 features: Tb_scaled, interaction, intercept
        X_model = np.column_stack([tb_scaled, interaction])
        model = LinearRegression().fit(X_model, y_s)
        
        A  = float(model.coef_[0])
        B  = float(model.coef_[1])
        C  = float(model.intercept_)

        logging.info(f"  {sensor}: A={A:.6f}  B={B:.6f}  C={C:.6f}")

        master_params[sensor] = {
            'Coeff_A_Tb':         A,
            'Coeff_B_TbTCWV': B,
            'Coeff_Intercept':    C,
            'Scaler_Tb_Center':   float(scaler.center_[0]),
            'Scaler_Tb_Scale':    float(scaler.scale_[0]),
            'Scaler_TCWV_Center': float(scaler.center_[1]),
            'Scaler_TCWV_Scale':  float(scaler.scale_[1]),
            'D_mean':             0.0,   # placeholder, filled after validation
            'D_std':              0.0,
        }

        # Apply back to training set for diagnostic plot
        full_X = training_df[feature_cols]
        valid_idx = full_X.dropna().index
        if not valid_idx.empty:
            X_v = scaler.transform(full_X.loc[valid_idx].values)
            tb_v = X_v[:, 0]
            tcwv_v = X_v[:, 1]
            interaction_v = tb_v * tcwv_v
            training_df.loc[valid_idx, f'retrieved_SST_{sensor}'] = (
                A * tb_v + B * interaction_v + C
            )

    params_df = pd.DataFrame.from_dict(master_params, orient='index')
    return params_df, scalers, training_df

# ============================================================================
# 2. Load and apply coefficients to in-situ validation data.
# ============================================================================

def validate_coefficients(validation_dir: Path,
                           params_df: pd.DataFrame,
                           scalers: dict) -> pd.DataFrame:
    """
    Apply dervied A, B, C coefficients to the in-situ validation matchups.
    Returns a DataFrame with a raw (pre-D) SST retrieval column added.
    """
    logging.info("--- PHASE 2: Applying coefficients to validation data ---")

    all_files = glob.glob(str(Path(validation_dir) / "CORA*.csv"))
    if not all_files:
        logging.error(f"No validation CSV files found in {validation_dir}.")
        return None

    df_list = []
    for f in all_files:
        try:
            df = pd.read_csv(f, low_memory=False, on_bad_lines='skip')
            # Normalise TCWV column name (called TQV in merra2)
            if "TQV" in df.columns and "TCWV" not in df.columns:
                df = df.rename(columns={"TQV": "TCWV"})

            # --- Quality filters ---
            df['TEMP'] = pd.to_numeric(df['TEMP'], errors='coerce')
            df['brightness_temp'] = pd.to_numeric(df['brightness_temp'], errors='coerce')
            df = df[(df['TEMP'] >= -1.8) & (df['TEMP'] <= 15.0)] # fitler to plausible SST range (-1.8 °C to 15 °C)
            df = df[(df['brightness_temp'] >= 271.35) & (df['brightness_temp'] <= 288.15)]
            if "ERA5_wind_speed" in df.columns:
                df['ERA5_wind_speed'] = pd.to_numeric(df['ERA5_wind_speed'], errors='coerce')
                df = df[df['ERA5_wind_speed'] >= 6.0]
            df = df[~df['TYPE'].isin(['FB', 'CT'])] # filter out ferrboxes and ctd files (should already be excluded in validation)
            if "pixel_count" in df.columns:
                df = df[df['pixel_count'].astype(float) > 5] # must have at least 5 pixels (100 m buffer at 30 m res = 9 pixels)

            # uncomment to filter out matchups with high within buffer temp variation
            # # stdev > 1.0 K (vs TIRS NEdT ~0.4 K) indicates spatial heterogeneity
            # # from contamination, not real SST variability
            # if "brightness_temp_stdev" in df.columns:
            #     df['brightness_temp_stdev'] = pd.to_numeric(
            #         df['brightness_temp_stdev'], errors='coerce'
            #     )
            #     df = df[df['brightness_temp_stdev'] < BT_STDEV_THRESHOLD]

            if not df.empty:
                df_list.append(df)
        except pd.errors.EmptyDataError:
            # Empty CSV — GEE batch had no Landsat matches, expected
            continue
        except Exception as e:
            logging.warning(f"Could not load {f}: {e}")
            continue

    if not df_list:
        logging.error("No valid validation data after loading all files.")
        return None

    val_df = pd.concat(df_list, ignore_index=True)
    val_df['sensor'] = val_df['LANDSAT_ID'].apply(get_sensor_from_id)
    val_df = val_df.dropna(subset=['sensor', 'TCWV', 'brightness_temp'])

    output_dfs = []
    for sensor in val_df['sensor'].unique():
        if sensor not in params_df.index:
            continue
        df_s = val_df[val_df['sensor'] == sensor].copy()
        if df_s.empty:
            continue

        row    = params_df.loc[sensor]
        scaler = scalers[sensor]

        X_val_scaled = scaler.transform(df_s[['brightness_temp', 'TCWV']].values)
        
        # Extract scaled terms
        tb_scaled = X_val_scaled[:, 0]
        tcwv_scaled = X_val_scaled[:, 1]
        interaction = tb_scaled * tcwv_scaled

        #  SST = A·Tb_scaled + B·(TCWV_scaled·Tb_scaled) + C (D derived later via k-fold validation)
        df_s['raw_SST_K'] = (
            row['Coeff_A_Tb']           * tb_scaled +
            row['Coeff_B_TbTCWV']  * interaction +
            row['Coeff_Intercept']
        )
        df_s['SST_in_situ_K'] = df_s['TEMP'] + 273.15
        output_dfs.append(df_s)

    if not output_dfs:
        return None

    final_df = pd.concat(output_dfs, ignore_index=True)

    if 'lon' not in final_df.columns or 'lat' not in final_df.columns:
        if 'LONGITUDE' in final_df.columns and 'LATITUDE' in final_df.columns:
            final_df['lon'] = pd.to_numeric(final_df['LONGITUDE'], errors='coerce')
            final_df['lat'] = pd.to_numeric(final_df['LATITUDE'],  errors='coerce')
        elif '.geo' in final_df.columns:
            def _extract_lonlat(geo_str):
                try:
                    coords = json.loads(geo_str).get('coordinates', [None, None])
                    return pd.Series({'lon': coords[0], 'lat': coords[1]})
                except Exception:
                    return pd.Series({'lon': None, 'lat': None})
            final_df = pd.concat([final_df, final_df['.geo'].apply(_extract_lonlat)], axis=1)

    return final_df

# ============================================================================
# 3. K-fold cross-validated derivation of D
# ============================================================================

def calculate_bias_kfold(val_df: pd.DataFrame):
    """
    Derive D via k-fold cross-validaton.

    Splits data into 5 folds per sensor. For each fold, derives D on the 
    other 4 folds with ±N_SIGMA outlier trimming, applies it to the fold that wasnt included, 
    and tags outliers in the test fold with the same ±N_SIGMA crieria. 
    Returns the combined (all 5 folds) test data with D applied, 
    all-data-with-D for background context in plots, per-sensor D values, and a stats table.
    """
    logging.info(f"--- PHASE 3: K-fold CV for D (K={CV_N_SPLITS}, ±{CV_N_SIGMA}σ trimming) ---")

    df_base = val_df.copy()
    df_base['Landsat']    = df_base['sensor']
    df_base['insitu_SST'] = df_base['SST_in_situ_K'] - 273.15
    df_base['raw_SST_C']  = df_base['raw_SST_K'] - 273.15

    test_fold_list  = []
    all_adj_list    = []
    optimal_D       = {}
    stats_records   = []

    kf = KFold(n_splits=CV_N_SPLITS, shuffle=True, random_state=CV_RANDOM_STATE)

    for sensor in sorted(df_base['Landsat'].unique()):
        df_s = df_base[df_base['Landsat'] == sensor].copy().reset_index(drop=True)
        n_total = len(df_s)

        if n_total < 2 * CV_N_SPLITS:
            logging.warning(
                f"  {sensor}: only {n_total} matchups; fewer than 2×K={2*CV_N_SPLITS}. "
                f"Using all data for D (no held-out test)."
            )
            error = df_s['raw_SST_C'] - df_s['insitu_SST']
            mu, sd = error.mean(), error.std()
            mask = (error >= mu - CV_N_SIGMA * sd) & (error <= mu + CV_N_SIGMA * sd)
            n_removed = (~mask).sum()
            if n_removed > 0:
                logging.info(f"  {sensor}: removed {n_removed} outliers from all data.")
            
            D = float(-(df_s[mask]['raw_SST_C'] - df_s[mask]['insitu_SST']).mean())
            df_s['adjusted_landsat_SST'] = df_s['raw_SST_C'] + D
            # Tag outliers on the (only) fold using the same ±N_SIGMA criri
            test_res = df_s['adjusted_landsat_SST'] - df_s['insitu_SST']
            mu_t, sd_t = test_res.mean(), test_res.std()
            df_s['is_outlier'] = (
                (test_res < mu_t - CV_N_SIGMA * sd_t) |
                (test_res > mu_t + CV_N_SIGMA * sd_t)
            )
            inliers_fb = df_s[~df_s['is_outlier']]
            optimal_D[sensor] = D
            test_fold_list.append(df_s)
            all_adj_list.append(df_s)

            raw_s = compute_validation_stats(inliers_fb['insitu_SST'].values, inliers_fb['raw_SST_C'].values)
            cor_s = compute_validation_stats(inliers_fb['insitu_SST'].values, inliers_fb['adjusted_landsat_SST'].values)
            stats_records.append({
                'Sensor': sensor, 'N': n_total, 'N_inliers': len(inliers_fb),
                'D_mean': D, 'D_std': np.nan, 'note': 'all-data-fallback',
                'raw_RMSE': raw_s['rmse'], 'raw_accuracy': raw_s['accuracy'],
                'raw_precision': raw_s['precision'], 'raw_bias': raw_s['bias'],
                'cor_RMSE': cor_s['rmse'], 'cor_accuracy': cor_s['accuracy'],
                'cor_precision': cor_s['precision'], 'cor_bias': cor_s['bias'],
            })
            logging.info(f"  {sensor}: D={D:.4f}°C (fallback, no held-out)")
            continue

        # --- K-fold loop ---
        D_per_fold     = []
        fold_test_dfs  = []

        for fold_idx, (cal_idx, test_idx) in enumerate(kf.split(df_s)):
            df_cal  = df_s.iloc[cal_idx].copy()
            df_test = df_s.iloc[test_idx].copy()

            # Derive D from calibration fold with ±1σ SYMMETRIC trimming
            error_cal = df_cal['raw_SST_C'] - df_cal['insitu_SST']
            mu, sd = error_cal.mean(), error_cal.std()
            cal_mask = (error_cal >= mu - CV_N_SIGMA * sd) & (error_cal <= mu + CV_N_SIGMA * sd)
            n_removed = (~cal_mask).sum()
            if n_removed > 0:
                logging.info(
                    f"  {sensor} fold {fold_idx+1}: removed {n_removed} outliers "
                    f"from calibration set ({100*n_removed/len(df_cal):.1f}%)."
                )

            # D corrects the raw retrieval to match in-situ:
            # raw_SST + D ≈ insitu_SST  =>  D = mean(insitu - raw) on cal fold
            D_fold = float((df_cal[cal_mask]['insitu_SST'] -
                            df_cal[cal_mask]['raw_SST_C']).mean())
            D_per_fold.append(D_fold)

            # Apply this fold's D to the held-out test fold
            df_test['adjusted_landsat_SST'] = df_test['raw_SST_C'] + D_fold
            df_test['fold'] = fold_idx
            fold_test_dfs.append(df_test)

        # --- Aggregate fold results ---
        D_mean = float(np.mean(D_per_fold))
        D_std  = float(np.std(D_per_fold, ddof=1))
        optimal_D[sensor] = D_mean

        df_test_all_sensor = pd.concat(fold_test_dfs, ignore_index=True)

        # Tag test-fold outliers with the same ±N_SIGMA criterion used on cal folds.
        test_residuals = (df_test_all_sensor['adjusted_landsat_SST']
                          - df_test_all_sensor['insitu_SST'])
        mu_t, sd_t = test_residuals.mean(), test_residuals.std()
        outlier_flag = (
            (test_residuals < mu_t - CV_N_SIGMA * sd_t) |
            (test_residuals > mu_t + CV_N_SIGMA * sd_t)
        )
        df_test_all_sensor['is_outlier'] = outlier_flag
        n_outliers = int(outlier_flag.sum())
        if n_outliers > 0:
            logging.info(
                f"  {sensor}: {n_outliers} test-fold outliers tagged "
                f"({100*n_outliers/len(df_test_all_sensor):.1f}%) — "
                f"greyed in plots, excluded from stats."
            )

        # Stats on inliers only 
        inliers = df_test_all_sensor[~outlier_flag]
        raw_s = compute_validation_stats(
            inliers['insitu_SST'].values,
            inliers['raw_SST_C'].values,
        )
        cor_s = compute_validation_stats(
            inliers['insitu_SST'].values,
            inliers['adjusted_landsat_SST'].values,
        )

        logging.info(
            f"  {sensor}: D={D_mean:.4f}±{D_std:.4f}°C  "
            f"Raw RMSE={raw_s['rmse']:.3f}°C  "
            f"Corrected RMSE={cor_s['rmse']:.3f}°C  "
            f"N_inliers={len(inliers)}"
        )

        stats_records.append({
            'Sensor': sensor, 'N': n_total, 'N_inliers': len(inliers),
            'D_mean': D_mean, 'D_std': D_std, 'note': f'{CV_N_SPLITS}-fold CV',
            'raw_RMSE': raw_s['rmse'], 'raw_accuracy': raw_s['accuracy'],
            'raw_precision': raw_s['precision'], 'raw_bias': raw_s['bias'],
            'cor_RMSE': cor_s['rmse'], 'cor_accuracy': cor_s['accuracy'],
            'cor_precision': cor_s['precision'], 'cor_bias': cor_s['bias'],
        })

        test_fold_list.append(df_test_all_sensor)

        # Apply mean D to ALL data for background scatter
        df_s_all = df_s.copy()
        df_s_all['adjusted_landsat_SST'] = df_s_all['raw_SST_C'] + D_mean
        all_adj_list.append(df_s_all)

    df_test = pd.concat(test_fold_list, ignore_index=True)
    df_all  = pd.concat(all_adj_list,  ignore_index=True)

    # KDE density for scatter plot 
    df_test['density'] = 0.0
    for sensor in df_test['Landsat'].unique():
        sensor_mask = df_test['Landsat'] == sensor
        inlier_mask = sensor_mask & ~df_test['is_outlier']
        subset = df_test[inlier_mask]
        if len(subset) < 5:
            continue
        try:
            vals   = np.vstack([subset['adjusted_landsat_SST'], subset['insitu_SST']])
            kernel = stats.gaussian_kde(vals)(vals)
            df_test.loc[inlier_mask, 'density'] = kernel
        except Exception:
            pass

    df_test = df_test.sort_values(by='density')

    stats_table = pd.DataFrame(stats_records).set_index('Sensor')
    return df_test, df_all, optimal_D, stats_table

# ============================================================================
# 6. Plotting — training performance
# ============================================================================

def plot_diagnostics(df: pd.DataFrame):
    cols = [c for c in ['sst_skin_merra2_k', 'sim_bt_L8', 'tcwv_merra2_kg_m2']
            if c in df.columns]
    if not cols:
        return
    plot_df = df[cols].dropna()
    fig, axes = plt.subplots(1, len(cols), figsize=(18, 5))
    for i, col in enumerate(cols):
        axes[i].hist(plot_df[col], bins=30, edgecolor='black', alpha=0.7)
        axes[i].set_title(col)
    plt.tight_layout()
    plt.savefig(DIAGNOSTIC_PLOT_FILE, dpi=150)
    plt.close(fig)
    logging.info(f"Saved diagnostics plot to {DIAGNOSTIC_PLOT_FILE}")

def plot_training_scatter(df: pd.DataFrame, output_path: Path):
    logging.info("--- Plotting Training Data Performance ---")
    sensors = [s for s in SENSORS_TO_PROCESS if f'sim_bt_{s}' in df.columns]
    if not sensors:
        return
    cols = 3
    rows = (len(sensors) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    axes = np.array(axes).flatten()
    for i, sensor in enumerate(sensors):
        ax = axes[i]
        col_ret = f'retrieved_SST_{sensor}'
        if col_ret not in df.columns:
            ax.axis('off')
            continue
        subset = df[['sst_skin_merra2_k', col_ret]].dropna()
        x = subset['sst_skin_merra2_k'] - 273.15
        y = subset[col_ret] - 273.15
        rmse = np.sqrt(mean_squared_error(x, y))
        bias = float(np.mean(y - x))
        ax.scatter(x, y, alpha=0.05, s=1, color='steelblue')
        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lims, lims, 'k--', alpha=0.75)
        ax.set_title(f"Training — Landsat {sensor[1:]}")
        ax.set_xlabel("MERRA-2 skin T (°C)")
        ax.set_ylabel("Retrieved SST (°C)")
        ax.text(0.05, 0.95,
                f"RMSE = {rmse:.3f}°C\nBias = {bias:.3f}°C\nN = {len(x):,}",
                transform=ax.transAxes, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    for j in range(len(sensors), len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    logging.info(f"Saved training scatter to {output_path}")

# ============================================================================
# 7. Plotting — validation scatter 
# ============================================================================

def _draw_scatter_with_marginals(ax, df_sensor, orig_data,
                                  density_min, density_max,
                                  panel_label, sensor_label,
                                  raw_rmse=None, cor_rmse=None,
                                  D_mean=None, region_label=None):
    # Split test-fold data into inliers and outliers
    if 'is_outlier' in df_sensor.columns:
        df_inliers  = df_sensor[~df_sensor['is_outlier']]
        df_outliers = df_sensor[ df_sensor['is_outlier']]
    else:
        df_inliers  = df_sensor
        df_outliers = df_sensor.iloc[:0]   # empty

    # Layer 1: all-data background (very light grey — context only)
    sns.scatterplot(data=orig_data, x='insitu_SST', y='adjusted_landsat_SST',
                    color='grey', alpha=0.15, ax=ax, legend=False)
    # Layer 2: test-fold outliers (medium grey — visible but not prominent)
    if not df_outliers.empty:
        ax.scatter(df_outliers['insitu_SST'], df_outliers['adjusted_landsat_SST'],
                   color='#888888', alpha=0.4, s=12, zorder=3, linewidths=0)
    # Layer 3: test-fold inliers (KDE density coloured)
    sns.scatterplot(data=df_inliers, x='insitu_SST', y='adjusted_landsat_SST',
                    hue='density', palette='YlGnBu_r', alpha=0.6, ax=ax, legend=False,
                    hue_norm=(density_min, density_max))

    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_ylabel('')
    ax.set_xlabel('')

    # Marginal histograms — inliers only
    ax_histx = ax.inset_axes([0, 1, 1, 0.15], sharex=ax)
    ax_histy = ax.inset_axes([1, 0, 0.15, 1], sharey=ax)

    ax_histx.hist(df_inliers['insitu_SST'], bins=30,
                  color='#1f78b4', alpha=0.4, edgecolor='k', linewidth=0.3)
    ax_histy.hist(df_inliers['adjusted_landsat_SST'], bins=30,
                  color='#1f78b4', alpha=0.4, orientation='horizontal',
                  edgecolor='k', linewidth=0.3)

    for _ax, orient in [(ax_histx, 'x'), (ax_histy, 'y')]:
        _ax.tick_params(axis='x', labelbottom=(orient == 'y'))
        _ax.tick_params(axis='y', labelleft=(orient == 'x'))
        _ax.spines['top'].set_visible(False)
        _ax.spines['right'].set_visible(orient == 'x')
        _ax.spines['left'].set_visible(False)
        _ax.spines['bottom'].set_visible(orient == 'y')

    ax_histx.yaxis.tick_right()
    ax_histx.yaxis.set_label_position('right')
    ax_histx.yaxis.set_tick_params(pad=-25)
    try:
        m = int(round(ax_histx.get_ylim()[1], -1))
        h = int(round(0.5 * m, -1))
        if m > 0:
            ax_histx.set_yticks([h, m])
            ax_histx.set_yticklabels([h, m])
    except Exception:
        pass

    ax_histy.tick_params(axis='x', rotation=-90)
    ax_histy.xaxis.set_tick_params(pad=-25)
    try:
        m = int(round(ax_histy.get_xlim()[1], -1))
        h = int(round(0.5 * m, -1))
        if m > 0:
            ax_histy.set_xticks([h, m])
            ax_histy.set_xticklabels([h, m])
    except Exception:
        pass

    # 1:1 line
    lims = [ax.get_xlim()[0], ax.get_xlim()[1]]
    ax.plot(lims, lims, 'k--', alpha=0.6, zorder=5)

    # Panel label and sensor name
    ax.text(0.05, 0.97, panel_label, transform=ax.transAxes,
            fontsize=14, va='top')
    ax.text(0.20, 0.97, sensor_label, transform=ax.transAxes,
            fontsize=14, va='top', fontweight='bold')

    # Stats annotation — inliers only
    if not df_inliers.empty:
        accuracy  = float(np.median(np.abs(df_inliers['adjusted_landsat_SST'] - df_inliers['insitu_SST'])))
        precision = float(np.median(np.abs(
            df_inliers['adjusted_landsat_SST'] - df_inliers['insitu_SST'] -
            np.median(df_inliers['adjusted_landsat_SST'] - df_inliers['insitu_SST'])
        )))
        n_in  = len(df_inliers)
        n_out = len(df_outliers)
    else:
        accuracy, precision, n_in, n_out = 0.0, 0.0, 0, 0

    lines = []
    if region_label is not None:
        lines.append(region_label)
    lines += [
        r'$\mathrm{Accuracy}=%.2f$°C' % accuracy,
        r'$\mathrm{Precision}=%.2f$°C' % precision,
    ]
    if cor_rmse is not None:
        lines.append(r'$\mathrm{RMSE}=%.2f$°C' % cor_rmse)
    if D_mean is not None:
        lines.append(r'$D=%.2f$°C' % D_mean)
    lines.append(r'$n=%d$' % n_in)
    if n_out > 0:
        lines.append(r'$n_\mathrm{out}=%d$' % n_out)

    ax.text(0.05, 0.87, '\n'.join(lines),
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

def draw_validation_scatter_grid(df_test, df_all, output_path, optimal_D,
                                  stats_table=None,
                                  nrows=2, ncols=2, figsize=(12, 10)):
    logging.info(f"--- Drawing validation scatter grid -> {output_path} ---")
    plt.rcParams['font.family'] = 'Arial'

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                              sharex=True, sharey=True)
    axes = np.array(axes).flatten()

    sensors = np.sort(df_test['Landsat'].unique())[::-1]

    density_min = df_test['density'].min() if 'density' in df_test.columns else 0
    density_max = df_test['density'].max() if 'density' in df_test.columns else 1

    for i, sensor in enumerate(sensors):
        if i >= len(axes):
            break
        ax = axes[i]
        df_s    = df_test[df_test['Landsat'] == sensor]
        df_orig = df_all[df_all['Landsat'] == sensor]

        raw_rmse = float(stats_table.loc[sensor, 'raw_RMSE']) if (stats_table is not None and sensor in stats_table.index) else None
        cor_rmse = float(stats_table.loc[sensor, 'cor_RMSE']) if (stats_table is not None and sensor in stats_table.index) else None
        D_mean   = optimal_D.get(sensor)

        _draw_scatter_with_marginals(
            ax, df_s, df_orig, density_min, density_max,
            panel_label=chr(97 + i),
            sensor_label=f"Landsat {sensor[1:]}",
            raw_rmse=raw_rmse, cor_rmse=cor_rmse, D_mean=D_mean,
            region_label='Arctic:',
        )

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    # Shared colourbar
    if 'density' in df_test.columns and not df_test.empty:
        norm = mcolors.BoundaryNorm(np.linspace(density_min, density_max, 10), ncolors=256)
        sm = plt.cm.ScalarMappable(cmap='YlGnBu_r', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.08)
        cbar.set_label('KDE Density', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.yaxis.set_major_formatter(OOMFormatter(-2, "%1.1f"))

    fig.text(0.5, 0.04, 'In-situ SST (°C)', ha='center', fontsize=16)
    fig.text(0.06, 0.5, 'Landsat SST (°C)', va='center',
             rotation='vertical', fontsize=16)

    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved validation scatter to {output_path}")

# ============================================================================
# 8. Plotting — Greenland subset scatter
# ============================================================================

def validation_plot_just_greenland(df_test, df_all, output_path, optimal_D,
                                    stats_table=None):
    logging.info("--- Plotting Greenland validation scatter ---")
    try:
        world = gpd.read_file(COUNTRY_SHP)
        world = world.rename(columns={'ADMIN': 'name'})
        greenland = world[world['name'] == 'Greenland'].copy()
    except Exception as e:
        logging.warning(f"Could not load Greenland shapefile: {e}. Skipping.")
        return

    def _to_gdf(df, crs_out="EPSG:3413"):
        df = df.dropna(subset=['lon', 'lat']).copy()
        gdf = gpd.GeoDataFrame(df,
                geometry=[Point(r['lon'], r['lat']) for _, r in df.iterrows()],
                crs="EPSG:4326")
        return gdf.to_crs(crs_out)

    target_crs = "EPSG:3413"
    greenland_proj = greenland.to_crs(target_crs)
    gdf_test = _to_gdf(df_test)
    gdf_all  = _to_gdf(df_all)

    gdf_test_gl = gpd.sjoin_nearest(gdf_test, greenland_proj, how='inner',
                                     max_distance=200_000, distance_col="dist")
    gdf_all_gl  = gpd.sjoin_nearest(gdf_all, greenland_proj, how='inner',
                                     max_distance=200_000, distance_col="dist")

    gdf_test_gl = gdf_test_gl.to_crs("EPSG:4326")
    gdf_all_gl  = gdf_all_gl.to_crs("EPSG:4326")

    # Re-compute density for the Greenland subset
    for sensor in gdf_test_gl['Landsat'].unique():
        mask = gdf_test_gl['Landsat'] == sensor
        subset = gdf_test_gl[mask]
        if len(subset) < 5:
            gdf_test_gl.loc[mask, 'density'] = 0.0
            continue
        try:
            vals = np.vstack([subset['adjusted_landsat_SST'], subset['insitu_SST']])
            gdf_test_gl.loc[mask, 'density'] = stats.gaussian_kde(vals)(vals)
        except Exception:
            gdf_test_gl.loc[mask, 'density'] = 0.0
    gdf_test_gl = gdf_test_gl.sort_values('density')

    plt.rcParams['font.family'] = 'Arial'
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    sensors   = np.sort(gdf_test_gl['Landsat'].unique())[::-1]

    density_min = gdf_test_gl['density'].min()
    density_max = gdf_test_gl['density'].max()

    for i, sensor in enumerate(sensors):
        if i >= 3:
            break
        ax = axes[i]
        ax.set_box_aspect(1)
        df_s    = gdf_test_gl[gdf_test_gl['Landsat'] == sensor]
        df_orig = gdf_all_gl[gdf_all_gl['Landsat'] == sensor]

        raw_rmse = float(stats_table.loc[sensor, 'raw_RMSE']) if (stats_table is not None and sensor in stats_table.index) else None
        cor_rmse = float(stats_table.loc[sensor, 'cor_RMSE']) if (stats_table is not None and sensor in stats_table.index) else None
        D_mean   = optimal_D.get(sensor)

        _draw_scatter_with_marginals(
            ax, df_s, df_orig, density_min, density_max,
            panel_label=chr(97 + i),
            sensor_label=f"Landsat {sensor[1:]}",
            raw_rmse=raw_rmse, cor_rmse=cor_rmse, D_mean=D_mean,
            region_label='Greenland:',
        )

    for j in range(len(sensors), 3):
        axes[j].axis('off')

    if not gdf_test_gl.empty:
        norm = mcolors.BoundaryNorm(np.linspace(density_min, density_max, 10), ncolors=256)
        sm = plt.cm.ScalarMappable(cmap='YlGnBu_r', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.05)
        cbar.set_label('KDE Density', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.yaxis.set_major_formatter(OOMFormatter(-2, "%1.1f"))

    fig.text(0.5, 0.02, 'In-situ SST (°C)', ha='center', fontsize=16)
    fig.text(0.08, 0.5, 'Landsat SST (°C)', va='center',
             rotation='vertical', fontsize=16)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved Greenland scatter to {output_path}")

# ============================================================================
# 9. Plotting — map and diagnostics
# ============================================================================

def plot_validation_subplots(df_test, df_all, output_path):
    logging.info(f"--- Plotting Map Diagnostics -> {output_path} ---")

    df = df_test.copy().dropna(subset=['lon', 'lat', 'adjusted_landsat_SST', 'insitu_SST'])
    original_df = df_all.copy().dropna(subset=['lon', 'lat'])

    try:
        vals = np.vstack([df['lon'], df['lat']])
        df['density'] = stats.gaussian_kde(vals)(vals)
        df = df.sort_values('density')
    except Exception:
        df['density'] = 1.0

    if 'TYPE' in df.columns:
        df['TYPE'] = df['TYPE'].replace({
            'DR': 'Drifter Buoy', 'CT': 'CTD', 'BT': 'Bottle',
            'MO': 'Mooring', 'PF': 'Profilers', 'TS': 'Shipboard TS',
        })
    else:
        df['TYPE'] = 'Unknown'

    obs_sensor_counts = df.groupby(['TYPE', 'Landsat']).size().unstack(fill_value=0)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 14
    fig = plt.figure(figsize=(15, 12))
    gs  = fig.add_gridspec(3, 2, height_ratios=[2, 1, 1])

    # Map
    ax0 = fig.add_subplot(gs[0:2, 0], projection=ccrs.NorthPolarStereo())
    ax0.set_extent([-180, 180, 59, 90], ccrs.PlateCarree())
    ax0.add_feature(cfeature.OCEAN, zorder=0, facecolor='lightblue', alpha=0.5)
    ax0.add_feature(cfeature.LAND,  zorder=0, edgecolor='k', facecolor='w')
    ax0.coastlines()
    ax0.scatter(original_df['lon'], original_df['lat'], color='k', s=70,
                transform=ccrs.PlateCarree(), alpha=0.1, label='All data')
    sc = ax0.scatter(df['lon'], df['lat'], c=df['density'], cmap='YlGnBu_r', s=60,
                     transform=ccrs.PlateCarree(), alpha=0.9,
                     edgecolor='k', linewidth=0.2, label='Filtered (test)')
    cb = fig.colorbar(sc, ax=ax0, orientation='vertical', fraction=0.046, pad=0.02)
    cb.formatter.set_powerlimits((0, 0))
    cb.set_label('Spatial Density (KDE)')
    ax0.text(0.80, 0.95, f'N: {len(df)}', transform=ax0.transAxes,
             va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax0.text(0.05, 0.95, '(a)', transform=ax0.transAxes, va='top')
    gl = ax0.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                       alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(60, 91, 10))

    # Bar chart — obs type vs sensor
    colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99']
    if obs_sensor_counts.shape[1] > len(colors):
        colors = list(plt.cm.tab10(np.linspace(0, 1, obs_sensor_counts.shape[1])))
    else:
        colors = colors[:obs_sensor_counts.shape[1]]
    cmap = mcolors.ListedColormap(colors)

    ax1 = fig.add_subplot(gs[1, 1])
    if not obs_sensor_counts.empty:
        obs_sensor_counts.plot(kind='bar', stacked=True, ax=ax1, edgecolor='black',
                               linewidth=0.5, colormap=cmap)
    ax1.set_ylabel('Frequency')
    ax1.legend(title='Landsat')
    ax1.text(0.05, 0.95, '(b)', transform=ax1.transAxes, va='top')
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_xlabel('Observation Type')

    ax2 = fig.add_subplot(gs[2, 1])
    df_years = df_test.copy().dropna(subset=['adjusted_landsat_SST', 'insitu_SST'])
    if 'TIME' in df_years.columns:
        df_years['Year'] = pd.to_datetime(df_years['TIME'], errors='coerce').dt.year
        pivot = df_years.pivot_table(index='Year', columns='Landsat', values='insitu_SST',
                                     aggfunc='count', fill_value=0)
        if not pivot.empty:
            pivot.index = pivot.index.astype(int)
            pivot.plot(kind='bar', stacked=True, ax=ax2, edgecolor='black',
                       linewidth=0.5, colormap=cmap)
    ax2.set_ylabel('Frequency')
    ax2.set_xlabel('')
    ax2.text(0.05, 0.95, '(c)', transform=ax2.transAxes, va='top')

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"Saved map diagnostics to {output_path}")

# ============================================================================
# 10. Plotting — residuals boxplot
# ============================================================================

def plot_residuals_boxplot(df_test, output_path):
    logging.info(f"--- Plotting Residuals Boxplot -> {output_path} ---")
    # Exclude ±1σ outliers — consistent with scatter plots and reported stats
    if 'is_outlier' in df_test.columns:
        df = df_test[~df_test['is_outlier']].copy()
    else:
        df = df_test.copy()
    df['residual'] = df['adjusted_landsat_SST'] - df['insitu_SST']
    df['SST_bin']  = pd.cut(df['adjusted_landsat_SST'],
                             bins=[-2, 0, 2, 4, 6, 8, 10, 12, 15],
                             labels=['-2–0','0–2','2–4','4–6',
                                     '6–8','8–10','10–12','12–15'])
    sensors   = np.sort(df['Landsat'].unique())[::-1]
    n_sensors = len(sensors)
    plt.rcParams['font.family'] = 'Arial'
    fig, axes = plt.subplots(1, n_sensors, figsize=(5 * n_sensors, 6), sharey=True)
    if n_sensors == 1:
        axes = [axes]
    for i, sensor in enumerate(sensors):
        ax      = axes[i]
        df_s    = df[df['Landsat'] == sensor].dropna(subset=['SST_bin', 'residual'])
        if df_s.empty:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_title(f'Landsat {sensor[1:]}')
            continue
        cats = df_s['SST_bin'].cat.categories
        data = [df_s[df_s['SST_bin'] == c]['residual'].values
                for c in cats if c in df_s['SST_bin'].values]
        labels = [str(c) for c in cats if c in df_s['SST_bin'].values]
        ax.boxplot(data, labels=labels, patch_artist=True,
                   boxprops=dict(facecolor='lightblue', alpha=0.7),
                   medianprops=dict(color='red', linewidth=2),
                   whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=1.5))
        ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_xlabel('Landsat SST bin (°C)', fontsize=12)
        if i == 0:
            ax.set_ylabel('Residual (Landsat − in-situ) (°C)', fontsize=12)
        ax.set_title(f'Landsat {sensor[1:]}', fontsize=14, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
        for j, d in enumerate(data):
            ax.text(j + 1, ax.get_ylim()[0], f'n={len(d)}',
                    ha='center', va='top', fontsize=8, alpha=0.7)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    logging.info(f"Saved residuals boxplot to {output_path}")

# ============================================================================
# 11. Main
# ============================================================================

def main():
    setup_logging()
    logging.info("=" * 70)
    logging.info("derive_SST-COEFFS-main.py")
    logging.info(f"K-fold CV for D: K={CV_N_SPLITS}, ±{CV_N_SIGMA}σ trimming, BT_stdev < {BT_STDEV_THRESHOLD} K")
    logging.info("=" * 70)

    # 1. Load simulation results
    try:
        simulation_df = pd.read_csv(SIMULATION_RESULTS_CSV)
        # Backwards-compatable rename for era5 → merra2 column names
        simulation_df = simulation_df.rename(columns={
            'sst_skin_era5_k':  'sst_skin_merra2_k',
            'tcwv_era5_kg_m2':  'tcwv_merra2_kg_m2',
        })
        simulation_df = simulation_df[simulation_df['sst_skin_merra2_k'] >= SST_MIN_KELVIN].copy()
        logging.info(f"Loaded {len(simulation_df):,} simulation samples.")
    except Exception as e:
        logging.error(f"Could not load simulation results: {e}")
        return

    # 2. Diagnostic plots
    plot_diagnostics(simulation_df)

    # 3. Phase 1: derive A, B, C
    params_df, scalers, training_df = derive_coefficients(simulation_df)
    if params_df.empty:
        logging.error("Coefficient derivation failed.")
        return

    # 4. Plot training performance
    # plot_training_scatter(training_df, TRAINING_PLOT_FILE)

    # 5. Phase 2: apply coefficients to validation matchups
    val_df = validate_coefficients(VALIDATION_DATA_DIR, params_df, scalers)
    if val_df is None or val_df.empty:
        logging.warning("No validation data available. Saving parameters with D=0.")
        params_df.to_csv(MASTER_PARAMS_CSV, index_label='Sensor')
        return

    # 6. Phase 3: k-fold derivation of D
    df_test, df_all, optimal_D, stats_table = calculate_bias_kfold(val_df)

    # 7. Print and save validation statistics
    logging.info("\n--- Validation Statistics Summary ---")
    with pd.option_context('display.float_format', '{:.4f}'.format):
        logging.info("\n" + stats_table.to_string())

    stats_table.to_csv(VALIDATION_STATS_CSV, index_label='Sensor')
    logging.info(f"Saved validation statistics to {VALIDATION_STATS_CSV}")

    # 8. Plots
    draw_validation_scatter_grid(
        df_test, df_all, VALIDATION_PLOT_FILE, optimal_D,
        stats_table=stats_table, nrows=2, ncols=2, figsize=(12, 10),
    )

    validation_plot_just_greenland(
        df_test, df_all, VALIDATION_PLOT_GREENLAND, optimal_D,
        stats_table=stats_table,
    )

    plot_validation_subplots(df_test, df_all, VALIDATION_MAP_FILE)

    plot_residuals_boxplot(df_test, VALIDATION_BOXPLOT_FILE)

    # 9. Update master params with cross-validated D
    for sensor, d_mean in optimal_D.items():
        if sensor in params_df.index:
            params_df.loc[sensor, 'D_mean'] = d_mean
            if sensor in stats_table.index:
                params_df.loc[sensor, 'D_std'] = float(stats_table.loc[sensor, 'D_std'])

    params_df.to_csv(MASTER_PARAMS_CSV, index_label='Sensor')
    logging.info(f"Master parameters saved to {MASTER_PARAMS_CSV}")

    logging.info("\n--- Final Algorithm Parameters ---")
    with pd.option_context('display.float_format', '{:.6f}'.format):
        logging.info("\n" + params_df.to_string())

    logging.info("\nDone.")

if __name__ == "__main__":
    main()
