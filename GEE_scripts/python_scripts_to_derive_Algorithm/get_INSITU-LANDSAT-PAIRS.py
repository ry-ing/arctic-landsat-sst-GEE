#!/usr/bin/env python
"""
get_INSITU-LANDSAT-PAIRS.py

Matches EasyCORA in-situ ocean observations with Landsat 4/5/7/8/9 scenes
via Google Earth Engine, exports matched pairs as CSV to Google Drive.

Inputs: EasyCORA NetCDF files from CMEMS, Landsat TOA/SR in GEE,
MERRA-2 and ERA5 for auxilary met variables.

Output: per-batch CSV files in Google Drive with in-situ + Landsat columns.

Key choices: ±3h match window, 100m buffer, depth ≤3m, FB/CT types excluded.
Restartable — progress tracked in processing_state.json.
"""

import argparse
import ee
import glob
import json
import logging
import os
import re
import signal
import sys
import tempfile
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
import xarray as xr
from tqdm import tqdm

# ============================================================================
# 1. Configuration
# ============================================================================

# *** SET THIS to the local path of your EasyCORA arctic NetCDF directory ***
# Download from CMEMS: cmems_obs-ins_glo_phy-temp-sal_my_easycora_irr
PARENT_DATA_DIR = Path(
    r'/path/to/cmems_obs-ins_glo_phy-temp-sal_my_easycora_irr/arctic/'
)

REPO_ROOT = Path(__file__).resolve().parent.parent
STATE_DIR = REPO_ROOT / 'data/validation_logs'
LOG_FILE  = REPO_ROOT / 'logs/gee_matchup.log'
STATE_FILE = STATE_DIR / 'processing_state.json'

print(f"State file: {STATE_FILE}")
print(f"Log file: {LOG_FILE}")

PREFILTER_CSV = None
# --- Matchup parameters ---
GDRIVE_FOLDER = 'Arctic_Landsat_SST-main/landsat-insitu-matches'
GEE_START_DATE = '1982-07-01'
GEE_END_DATE = datetime.today().strftime('%Y-%m-%d')   # covers full EasyCORA archive
INSITU_START_DATE = pd.Timestamp(GEE_START_DATE)
INSITU_END_DATE = pd.Timestamp(GEE_END_DATE)
 
MATCH_WINDOW_HOURS = 3
GEOMETRY_BUFFER_METERS = 100
 
# --- In-situ QC ---
SURFACE_DEPTH_MAX = 3.0
EXCLUDED_TYPES: Set[str] = {'FB', 'CT'}
 
# --- Batching (v3) ---
BATCH_FILE_LIMIT = 50            # accumulate up to 50 files per GEE task
MAX_POINTS_PER_BATCH = 2000      # safety cap on total points per batch
 
# --- Task queue throttling ---
TASK_QUEUE_SOFT_LIMIT = 2500     # sleep if pending tasks above this
TASK_QUEUE_CHECK_EVERY = 10      # check queue every N batches
TASK_QUEUE_SLEEP_SECONDS = 60
 
# --- Retry 
FILE_READ_MAX_RETRIES = 3
FILE_READ_BACKOFF_SECONDS = 5
GEE_TASK_SUBMIT_RETRIES = 2
 
# --- Export 
EXPORT_COLUMNS = [
    'TIME', 'LATITUDE', 'LONGITUDE', 'TYPE', 'SOURCE_FILE',
    'DEPH', 'DEPH_min', 'DEPH_max', 'n_depths',
    'TEMP', 'TEMP_std',
    'LANDSAT_ID', 'brightness_temp', 'brightness_temp_stdev',
    'pixel_count', 'time_diff_seconds',
    'TQV', 'ERA5_temp_2m', 'ERA5_skin_temp', 'ERA5_wind_speed',
]
 
# ============================================================================
# 2. GEE Collection definitions
# ============================================================================
 
COLLECTION = {
    'L4': {
        'TOA': 'LANDSAT/LT04/C02/T1_TOA',
        'SR':  'LANDSAT/LT04/C02/T1_L2',
        'TIR': ['B6'],
        'VISW_SR':  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
        'VISW_TOA': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'],
        'has_pan': False, 'has_cirrus': False,
    },
    'L5': {
        'TOA': 'LANDSAT/LT05/C02/T1_TOA',
        'SR':  'LANDSAT/LT05/C02/T1_L2',
        'TIR': ['B6'],
        'VISW_SR':  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
        'VISW_TOA': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'],
        'has_pan': False, 'has_cirrus': False,
    },
    'L7': {
        'TOA': 'LANDSAT/LE07/C02/T1_TOA',
        'SR':  'LANDSAT/LE07/C02/T1_L2',
        'TIR': ['B6_VCID_1', 'B6_VCID_2'],
        'VISW_SR':  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
        'VISW_TOA': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B8'],
        'has_pan': True, 'has_cirrus': False,
    },
    'L8': {
        'TOA': 'LANDSAT/LC08/C02/T1_TOA',
        'SR':  'LANDSAT/LC08/C02/T1_L2',
        'TIR': ['B10', 'B11'],
        'VISW_SR':  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
        'VISW_TOA': ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9'],
        'has_pan': True, 'has_cirrus': True,
    },
    'L9': {
        'TOA': 'LANDSAT/LC09/C02/T1_TOA',
        'SR':  'LANDSAT/LC09/C02/T1_L2',
        'TIR': ['B10', 'B11'],
        'VISW_SR':  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
        'VISW_TOA': ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9'],
        'has_pan': True, 'has_cirrus': True,
    },
}
 
# ============================================================================
# 3. shutdown
# ============================================================================
 
STOP_REQUESTED = False
 
def _handle_signal(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True
    logging.warning(f"Received signal {signum}; finishing current batch then stopping.")
 
signal.signal(signal.SIGTERM, _handle_signal)
signal.signal(signal.SIGINT, _handle_signal)
 
# ============================================================================
# 4. State management
# ============================================================================
#
# Status values:
#   'exported'       -> done, data confirmed in Drive
#   'empty'          -> no surface obs, or no Landsat match found
#   'skipped_type'   -> excluded instrumnet type (FB/CT)
#   'submitted'      -> queued in GEE, not yet verified
#   'read_failed'    -> transient; retried on restart
#   'errored'        -> transient; retried on restart

TERMINAL_STATUSES = {'exported', 'empty', 'skipped_type'}
# 'submitted' skipped by main loop (already in flight) but NOT terminal —
# the --verify pass resolves it to exported/empty after Drive syncs.
SKIP_STATUSES = TERMINAL_STATUSES | {'submitted'}
 
def load_state() -> Dict[str, dict]:
    if not STATE_FILE.exists():
        return {}
    try:
        with open(STATE_FILE) as f:
            return json.load(f)
    except Exception as e:
        logging.error(f"Could not read state file {STATE_FILE}: {e}. Starting fresh.")
        backup = STATE_FILE.with_suffix(f'.corrupt.{int(time.time())}.json')
        try:
            STATE_FILE.rename(backup)
            logging.error(f"Corrupt state moved to {backup}.")
        except Exception:
            pass
        return {}
 
def save_state_atomic(state: Dict[str, dict]) -> None:
    STATE_FILE.parent.mkdir(exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(
        prefix='state_', suffix='.json.tmp', dir=str(STATE_FILE.parent)
    )
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(state, f, indent=2, sort_keys=True)
        os.replace(tmp_path, STATE_FILE)
    except Exception as e:
        logging.error(f"Failed atomic state write: {e}")
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
 
def mark_file(state: Dict[str, dict], filepath: Path, status: str, **kwargs) -> None:
    record = {
        'status':    status,
        'timestamp': datetime.utcnow().isoformat(timespec='seconds'),
    }
    record.update(kwargs)
    state[str(filepath.resolve())] = record
 
def mark_many_and_save(state: Dict[str, dict], files: List[Path], status: str, **kwargs) -> None:
    for fp in files:
        mark_file(state, fp, status=status, **kwargs)
    save_state_atomic(state)
 
def is_terminal(state: Dict[str, dict], filepath: Path) -> bool:
    """True if file should be skipped by main loop. Includes 'submitted' — use --verify to resolve."""
    rec = state.get(str(filepath.resolve()))
    if rec is None:
        return False
    return rec.get('status') in SKIP_STATUSES
 
# ============================================================================
# 5. GEE masking / loading
# ============================================================================
 
def _set_time_diff(target_date: ee.Date):
    def wrapper(img: ee.Image) -> ee.Image:
        return img.set(
            'time_diff',
            ee.Number(img.get('system:time_start'))
              .subtract(target_date.millis()).abs()
        )
    return wrapper
 
 
def mask_water(image: ee.Image, landsat: str) -> ee.Image:
    qa = image.select('QA_PIXEL')
 
    is_water = qa.bitwiseAnd(1 << 7).neq(0)
 
    bad = (qa.bitwiseAnd(1 << 0).eq(0)                    # fill
             .And(qa.bitwiseAnd(1 << 1).eq(0))            # dilated cloud
             .And(qa.bitwiseAnd(1 << 3).eq(0))            # cloud
             .And(qa.bitwiseAnd(1 << 4).eq(0))            # cloud shadow
             .And(qa.bitwiseAnd(1 << 5).eq(0)))           # snow
 
    cloud_conf  = qa.rightShift(8).bitwiseAnd(3)
    shadow_conf = qa.rightShift(10).bitwiseAnd(3)
    snow_conf   = qa.rightShift(12).bitwiseAnd(3)
    bad = bad.And(cloud_conf.lt(2)).And(shadow_conf.lt(2)).And(snow_conf.lt(2))
 
    if COLLECTION[landsat]['has_cirrus']:
        bad = bad.And(qa.bitwiseAnd(1 << 2).eq(0))
        cirrus_conf = qa.rightShift(14).bitwiseAnd(3)
        bad = bad.And(cirrus_conf.lt(2))
 
    combined_mask = is_water.And(bad)
 
    bt_masked = image.select('brightness_temp').updateMask(combined_mask)
    return image.addBands(bt_masked.rename('brightness_temp_masked'), overwrite=True)
 
def get_landsat_data(landsat: str, date_start: str, date_end: str) -> ee.ImageCollection:
    cfg = COLLECTION[landsat]
    tir_primary = cfg['TIR'][0]
 
    toa = (ee.ImageCollection(cfg['TOA'])
             .filter(ee.Filter.date(date_start, date_end))
             .map(lambda img: (img.select(cfg['VISW_TOA'])
                                  .addBands(img.select([tir_primary], ['brightness_temp'])))))
 
    sr = (ee.ImageCollection(cfg['SR'])
            .filter(ee.Filter.date(date_start, date_end))
            .select(cfg['VISW_SR']))
 
    joined = sr.combine(toa, True)
    return joined.map(lambda img: mask_water(img, landsat))
 
def load_landsat_master(date_start: str, date_end: str) -> ee.ImageCollection:
    logging.info("Building Landsat master collection...")
    colls = [get_landsat_data(s, date_start, date_end)
             for s in ['L9', 'L8', 'L7', 'L5', 'L4']]
    master = colls[0]
    for c in colls[1:]:
        master = master.merge(c)
    return master.sort('system:time_start')
 
def load_merra2_collection(date_start: str, date_end: str) -> ee.ImageCollection:
    return (ee.ImageCollection("NASA/GSFC/MERRA/slv/2")
              .filter(ee.Filter.date(date_start, date_end))
              .select('TQV'))
 
def load_era5_collection(date_start: str, date_end: str) -> ee.ImageCollection:
    return (ee.ImageCollection("ECMWF/ERA5/HOURLY")
              .filter(ee.Filter.date(date_start, date_end))
              .select(['temperature_2m', 'skin_temperature',
                       'u_component_of_wind_10m', 'v_component_of_wind_10m']))
 
# ============================================================================
# 6. Matchup and extraction (server-side)
# ============================================================================
 
def flag_matching_obs(feature: ee.Feature, collection: ee.ImageCollection) -> ee.Feature:
    date = ee.Date(feature.get('TIME'))
    sdate = date.advance(-MATCH_WINDOW_HOURS, 'hour')
    edate = date.advance( MATCH_WINDOW_HOURS, 'hour')
    location = feature.geometry().buffer(GEOMETRY_BUFFER_METERS)
 
    image = (collection
             .filterBounds(location)
             .filterDate(sdate, edate)
             .map(_set_time_diff(date))
             .sort('time_diff')
             .first())
 
    exists = ee.Algorithms.If(image, True, False)
    return feature.set('LANDSAT_exists', exists)
 
def extract_matchup(feature: ee.Feature,
                    landsat_coll: ee.ImageCollection,
                    merra_coll:   ee.ImageCollection,
                    era5_coll:    ee.ImageCollection) -> ee.Feature:
    date = ee.Date(feature.get('TIME'))
    sdate = date.advance(-MATCH_WINDOW_HOURS, 'hour')
    edate = date.advance( MATCH_WINDOW_HOURS, 'hour')
 
    landsat_area = feature.geometry().buffer(GEOMETRY_BUFFER_METERS)
    pt = feature.geometry()
    time_differ = _set_time_diff(date)
 
    landsat_image = (landsat_coll
                     .filterBounds(landsat_area)
                     .filterDate(sdate, edate)
                     .map(time_differ)
                     .sort('time_diff')
                     .first())
 
    time_diff_s = ee.Number(landsat_image.get('time_diff')).divide(1000)
 
    bt_stats = (landsat_image.select('brightness_temp_masked')
                .reduceRegion(
                    reducer=(ee.Reducer.mean()
                             .combine(ee.Reducer.count(),  '', True)
                             .combine(ee.Reducer.stdDev(), '', True)),
                    geometry=landsat_area,
                    scale=100,
                    maxPixels=1e6,
                ))
 
    brightness_temp       = bt_stats.get('brightness_temp_masked_mean')
    brightness_temp_stdev = bt_stats.get('brightness_temp_masked_stdDev')
    pixel_count           = bt_stats.get('brightness_temp_masked_count')
    landsat_id            = landsat_image.get('LANDSAT_PRODUCT_ID')
 
    merra_image = (merra_coll
                   .filterBounds(pt)
                   .filterDate(sdate, edate)
                   .map(time_differ)
                   .sort('time_diff')
                   .first())
    tqv = merra_image.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=pt, scale=50000
    ).get('TQV')
 
    era5_image = (era5_coll
                  .filterBounds(pt)
                  .filterDate(sdate, edate)
                  .map(time_differ)
                  .sort('time_diff')
                  .first())
    era5_stats = era5_image.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=pt, scale=25000
    )
    t2m    = era5_stats.get('temperature_2m')
    tskin  = era5_stats.get('skin_temperature')
    u_wind = era5_stats.get('u_component_of_wind_10m')
    v_wind = era5_stats.get('v_component_of_wind_10m')
    wind_speed = ee.Number(u_wind).pow(2).add(ee.Number(v_wind).pow(2)).sqrt()
 
    return feature.set(
        'brightness_temp',        brightness_temp,
        'brightness_temp_stdev',  brightness_temp_stdev,
        'pixel_count',            pixel_count,
        'LANDSAT_ID',             landsat_id,
        'time_diff_seconds',      time_diff_s,
        'TQV',                    tqv,
        'ERA5_temp_2m',           t2m,
        'ERA5_skin_temp',         tskin,
        'ERA5_wind_speed',        wind_speed,
    )
 
# ============================================================================
# 7. In-situ loading
# ============================================================================
 
def _read_netcdf_with_retry(path: Path) -> Optional[xr.Dataset]:
    last_err = None
    for attempt in range(FILE_READ_MAX_RETRIES):
        try:
            return xr.open_dataset(path, engine='netcdf4')
        except (OSError, RuntimeError) as e:
            last_err = e
            wait = FILE_READ_BACKOFF_SECONDS * (2 ** attempt)
            logging.warning(
                f"Read error on {path.name} (attempt {attempt+1}/{FILE_READ_MAX_RETRIES}): "
                f"{e!r}; sleeping {wait}s."
            )
            time.sleep(wait)
        except Exception as e:
            logging.warning(f"Non-retryable error opening {path.name}: {e!r}")
            return None
    logging.error(f"Giving up on {path.name} after {FILE_READ_MAX_RETRIES} retries. Last error: {last_err!r}")
    return None
 
def load_insitu_data(path: Path) -> Optional[pd.DataFrame]:
    ds = _read_netcdf_with_retry(path)
    if ds is None:
        return None
 
    try:
        try:
            t0 = pd.Timestamp(ds.TIME.values[0])
        except Exception:
            t0 = pd.Timestamp(pd.to_datetime(ds.TIME.values)[0])
        if not (INSITU_START_DATE <= t0 <= INSITU_END_DATE):
            return pd.DataFrame()
 
        m = re.search(r'_(\w{2})\.nc$', path.name)
        itype = m.group(1) if m else 'N/A'
 
        if itype in EXCLUDED_TYPES:
            return pd.DataFrame()
 
        if 'DEPH' not in ds:
            if 'PRES' in ds and 'DENS' in ds:
                ds['DEPH'] = (ds['PRES'] * 1000.0) / (ds['DENS'] * 9.81)
            else:
                return None
 
        df = ds[['DEPH', 'TEMP', 'TIME', 'LATITUDE', 'LONGITUDE']].to_dataframe()
        df['TYPE'] = itype
    finally:
        try:
            ds.close()
        except Exception:
            pass
 
    df = df.dropna(subset=['TEMP', 'DEPH', 'LATITUDE', 'LONGITUDE', 'TIME'])
    df = df[df['DEPH'] <= SURFACE_DEPTH_MAX]
    if df.empty:
        return df
 
    df['TIME'] = pd.to_datetime(df['TIME'])
    df = df[(df['LATITUDE'].between(-90, 90)) & (df['LONGITUDE'].between(-180, 180))]
    df = df[~((df['LATITUDE'] == 0) & (df['LONGITUDE'] == 0))]
    if df.empty:
        return df
 
    agg = (df.groupby(['TIME', 'LATITUDE', 'LONGITUDE', 'TYPE'], as_index=False)
             .agg(TEMP=('TEMP', 'mean'),
                  TEMP_std=('TEMP', 'std'),
                  DEPH=('DEPH', 'mean'),
                  DEPH_min=('DEPH', 'min'),
                  DEPH_max=('DEPH', 'max'),
                  n_depths=('TEMP', 'size')))
    agg['TEMP_std'] = agg['TEMP_std'].fillna(0.0)
    return agg
 
def pandas_to_ee(df: pd.DataFrame) -> ee.FeatureCollection:
    features = []
    for _, row in df.iterrows():
        props = {
            'TIME':         row['TIME'].strftime('%Y-%m-%dT%H:%M:%S'),
            'DEPH':         float(row['DEPH']),
            'DEPH_min':     float(row['DEPH_min']),
            'DEPH_max':     float(row['DEPH_max']),
            'TEMP':         float(row['TEMP']),
            'TEMP_std':     float(row.get('TEMP_std', 0.0)),
            'n_depths':     int(row.get('n_depths', 1)),
            'TYPE':         str(row.get('TYPE', 'N/A')),
            'SOURCE_FILE':  str(row.get('source_file', '')),
            'LATITUDE':     float(row['LATITUDE']),
            'LONGITUDE':    float(row['LONGITUDE']),
        }
        geom = ee.Geometry.Point([float(row['LONGITUDE']), float(row['LATITUDE'])])
        features.append(ee.Feature(geom, props))
    return ee.FeatureCollection(features)
 
# ============================================================================
# 8. Task submission and queue throttling
# ============================================================================
 
def count_pending_tasks() -> int:
    """Count pending/running tasks; tries newer listOperations API first."""
    try:
        ops = ee.data.listOperations()
        return sum(
            1 for op in ops
            if op.get('metadata', {}).get('state') in ('PENDING', 'RUNNING')
        )
    except AttributeError:
        pass
    except Exception as e:
        logging.warning(f"listOperations() failed: {e!r}")
    # fallback for older GEE API
    try:
        tasks = ee.data.getTaskList()
        return sum(1 for t in tasks if t.get('state') in ('READY', 'RUNNING'))
    except Exception as e:
        logging.warning(f"Could not query task queue: {e!r}; assuming 0.")
        return 0
 
def throttle_if_needed() -> None:
    while not STOP_REQUESTED:
        pending = count_pending_tasks()
        if pending <= TASK_QUEUE_SOFT_LIMIT:
            return
        logging.info(
            f"Task queue has {pending} pending (soft limit {TASK_QUEUE_SOFT_LIMIT}); "
            f"sleeping {TASK_QUEUE_SLEEP_SECONDS}s."
        )
        time.sleep(TASK_QUEUE_SLEEP_SECONDS)
 
def submit_batch(df: pd.DataFrame,
                 landsat_coll: ee.ImageCollection,
                 merra_coll:   ee.ImageCollection,
                 era5_coll:    ee.ImageCollection,
                 batch_id:     str) -> Optional[str]:
    """Submit one batch containg points from many files. Returns task ID or None on failure."""
    if df.empty:
        return None
 
    fc = pandas_to_ee(df)
    flagged = fc.map(lambda f: flag_matching_obs(f, landsat_coll))
    matched = flagged.filter(ee.Filter.eq('LANDSAT_exists', True))
    extracted = matched.map(
        lambda f: extract_matchup(f, landsat_coll, merra_coll, era5_coll)
    )
    final = extracted.filter(ee.Filter.neq('brightness_temp', None))
 
    description = f"CORA_Landsat_Validation_{batch_id}"
    last_err = None
    for attempt in range(GEE_TASK_SUBMIT_RETRIES + 1):
        try:
            task = ee.batch.Export.table.toDrive(
                collection=final,
                description=description,
                folder=GDRIVE_FOLDER,
                fileFormat='CSV',
                selectors=EXPORT_COLUMNS,
            )
            task.start()
            logging.info(f"Batch {batch_id}: submitted task {task.id} "
                         f"({len(df)} points).")
            return task.id
        except Exception as e:
            last_err = e
            wait = 2 ** attempt
            logging.warning(f"Batch {batch_id}: task.start() failed ({e!r}); retry in {wait}s.")
            time.sleep(wait)
    logging.error(f"Batch {batch_id}: gave up submitting. Last error: {last_err!r}")
    return None
 
# ============================================================================
# 9. Verification pass
# ============================================================================
 
def get_all_task_statuses() -> Dict[str, dict]:
    try:
        ops = ee.data.listOperations()
        out = {}
        for op in ops:
            name = op.get('name', '')
            tid = name.split('/')[-1] if '/' in name else name
            meta = op.get('metadata', {})
            out[tid] = {
                'state': meta.get('state', 'UNKNOWN'),
                'description': meta.get('description', ''),
                'destination_uris': meta.get('destinationUris', []),
            }
        return out
    except AttributeError:
        pass
    except Exception as e:
        logging.warning(f"listOperations() failed: {e}")
    try:
        tasks = ee.data.getTaskList()
        return {t['id']: t for t in tasks}
    except Exception as e:
        logging.error(f"Could not list tasks: {e}")
        return {}
 
def read_csv_source_files_from_drive(csv_dir: Path) -> Set[str]:
    """Read all CSVs and collect the set of SOURCE_FILE values present."""
    if not csv_dir.exists():
        logging.error(f"CSV dir does not exist: {csv_dir}")
        return set()
    seen = set()
    csvs = list(csv_dir.glob("CORA_Landsat_Validation_*.csv"))
    logging.info(f"Reading {len(csvs)} CSVs from {csv_dir}.")
    for csv_path in tqdm(csvs, desc='Reading CSVs'):
        try:
            df = pd.read_csv(csv_path, usecols=['SOURCE_FILE'])
            seen.update(df['SOURCE_FILE'].astype(str).unique())
        except ValueError:
            # CSV may be empty (no matches in batch)
            continue
        except Exception as e:
            logging.warning(f"Could not read {csv_path.name}: {e}")
    return seen
 
def verify(csv_dir: Optional[Path]) -> None:
    if csv_dir is None:
        logging.error("--verify requires --csv-dir pointing to a local copy "
                      "of the Drive folder.")
        return
 
    state = load_state()
    submitted = [
        Path(p) for p, rec in state.items()
        if rec.get('status') == 'submitted'
    ]
    if not submitted:
        logging.info("No 'submitted' files to verify.")
        return
 
    logging.info(f"Verifying {len(submitted)} submitted files against CSVs in {csv_dir}.")
 
    source_files_present = read_csv_source_files_from_drive(csv_dir)
    logging.info(
        f"Found {len(source_files_present)} distinct SOURCE_FILE values in CSVs."
    )
 
    task_statuses = get_all_task_statuses()
 
    n_exp, n_emp, n_run, n_fail = 0, 0, 0, 0
    for fp in submitted:
        rec = state[str(fp.resolve())]
        tid = rec.get('task_id')
        task_state = task_statuses.get(tid, {}).get('state', 'UNKNOWN')
 
        if task_state in ('READY', 'RUNNING', 'PENDING'):
            n_run += 1
            continue
 
        if task_state == 'FAILED':
            mark_file(state, fp, status='errored',
                      error='GEE task FAILED',
                      previous_batch_id=rec.get('batch_id'),
                      previous_task_id=tid)
            n_fail += 1
            continue
 
        # COMPLETED or UNKNOWN: rely on CSV presence
        if fp.name in source_files_present:
            mark_file(state, fp, status='exported',
                      batch_id=rec.get('batch_id'),
                      task_id=tid)
            n_exp += 1
        else:
            mark_file(state, fp, status='empty',
                      batch_id=rec.get('batch_id'),
                      task_id=tid,
                      note='no matched rows in exported CSV')
            n_emp += 1
 
    save_state_atomic(state)
    logging.info(
        f"Verification: exported={n_exp}, empty={n_emp}, "
        f"still_running={n_run}, failed={n_fail}"
    )
 
# ============================================================================
# 10. Main processing
# ============================================================================
 
def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler(sys.stdout),
        ],
    )
 
def get_file_list() -> List[Path]:
    logging.info(f"Globbing .nc files in {PARENT_DATA_DIR} ...")
    all_files = [Path(f) for f in glob.glob(str(PARENT_DATA_DIR / '**/*.nc'), recursive=True)]
    logging.info(f"Found {len(all_files)} .nc files.")
 
    if PREFILTER_CSV and PREFILTER_CSV.exists():
        logging.info(f"Applying pre-filter from {PREFILTER_CSV} ...")
        try:
            prefilter = pd.read_csv(PREFILTER_CSV)
            keep_stems = set(prefilter['insitu_filename'].astype(str))
            all_files = [f for f in all_files if f.stem in keep_stems]
            logging.info(f"After pre-filter: {len(all_files)} files.")
        except Exception as e:
            logging.warning(f"Pre-filter load failed ({e!r}); keeping all files.")
 
    return sorted(all_files)
 
def main_process() -> None:
    logging.info("=" * 72)
    logging.info("Starting get_INSITU-LANDSAT-PAIRS.py (v3 - batched)")
    logging.info(f"Match window: +/- {MATCH_WINDOW_HOURS} hours")
    logging.info(f"Batch file limit: {BATCH_FILE_LIMIT}, max points/batch: {MAX_POINTS_PER_BATCH}")
    logging.info(f"Task queue soft limit: {TASK_QUEUE_SOFT_LIMIT}")
    logging.info(f"State file: {STATE_FILE}")
    logging.info("=" * 72)
 
    try:
        ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
    except Exception as e:
        logging.error(f"Earth Engine auth failed: {e}")
        logging.error("Run 'earthengine authenticate' first.")
        sys.exit(1)
 
    try:
        landsat_master = load_landsat_master(GEE_START_DATE, GEE_END_DATE)
        merra2_tqv     = load_merra2_collection(GEE_START_DATE, GEE_END_DATE)
        era5           = load_era5_collection(GEE_START_DATE, GEE_END_DATE)
    except Exception as e:
        logging.error(f"Failed to build GEE collections: {e}")
        sys.exit(1)
 
    state = load_state()
 
    status_counts: Dict[str, int] = {}
    for rec in state.values():
        status_counts[rec.get('status', '?')] = status_counts.get(rec.get('status', '?'), 0) + 1
    logging.info(
        "Existing state: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items()))
    )
 
    all_files = get_file_list()
    pending = [f for f in all_files if not is_terminal(state, f)]
    logging.info(f"Files pending: {len(pending)}/{len(all_files)}")
 
    if not pending:
        logging.info("Nothing to do.")
        return
 
    run_tag = datetime.utcnow().strftime('%Y%m%d-%H%M%S')
    batch_counter = 0
    batches_since_queue_check = 0
 
    acc_dfs: List[pd.DataFrame] = []
    acc_files: List[Path] = []
    acc_points = 0
 
    def flush_batch() -> None:
        nonlocal batch_counter, acc_dfs, acc_files, acc_points
        nonlocal batches_since_queue_check
 
        if not acc_dfs:
            return
 
        batches_since_queue_check += 1
        if batches_since_queue_check >= TASK_QUEUE_CHECK_EVERY:
            throttle_if_needed()
            batches_since_queue_check = 0
 
        batch_counter += 1
        batch_id = f"run{run_tag}-{batch_counter:04d}"
        combined = pd.concat(acc_dfs, ignore_index=True)
 
        logging.info(
            f"Flushing batch {batch_id}: {len(acc_files)} files, {len(combined)} points."
        )
        task_id = None
        try:
            task_id = submit_batch(
                combined, landsat_master, merra2_tqv, era5, batch_id=batch_id
            )
        except Exception as e:
            logging.error(f"Batch {batch_id} submission crashed: {e!r}")
            logging.error(traceback.format_exc())
 
        if task_id:
            mark_many_and_save(
                state, acc_files, status='submitted',
                batch_id=batch_id, task_id=task_id,
                n_points_in_batch=int(len(combined)),
                n_files_in_batch=len(acc_files),
            )
        else:
            mark_many_and_save(
                state, acc_files, status='errored',
                error='batch submission failed',
                batch_id=batch_id,
            )
 
        acc_dfs, acc_files, acc_points = [], [], 0
 
    for file_path in tqdm(pending, desc='Processing files'):
        if STOP_REQUESTED:
            logging.warning("Stop requested; flushing pending batch then exiting.")
            break
 
        try:
            insitu_df = load_insitu_data(file_path)
        except Exception as e:
            logging.error(f"{file_path.name}: unexpected load error: {e!r}")
            mark_file(state, file_path, status='errored',
                      error=f"load crash: {e!r}"[:500])
            save_state_atomic(state)
            continue
 
        if insitu_df is None:
            mark_file(state, file_path, status='read_failed')
            save_state_atomic(state)
            continue
 
        if insitu_df.empty:
            m = re.search(r'_(\w{2})\.nc$', file_path.name)
            itype = m.group(1) if m else 'N/A'
            status = 'skipped_type' if itype in EXCLUDED_TYPES else 'empty'
            mark_file(state, file_path, status=status, instrument_type=itype)
            save_state_atomic(state)
            continue
 
        insitu_df['source_file'] = file_path.name
 
        acc_dfs.append(insitu_df)
        acc_files.append(file_path)
        acc_points += len(insitu_df)
 
        if (len(acc_files) >= BATCH_FILE_LIMIT
                or acc_points >= MAX_POINTS_PER_BATCH):
            flush_batch()
 
    if acc_files:
        logging.info("Flushing final partial batch.")
        flush_batch()
 
    status_counts = {}
    for rec in state.values():
        status_counts[rec.get('status', '?')] = status_counts.get(rec.get('status', '?'), 0) + 1
    logging.info(
        "Run summary: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items()))
    )
    logging.info(
        "After GEE tasks complete and CSVs sync locally, run with "
        "`--verify --csv-dir /path/to/csvs` to resolve 'submitted' entries."
    )
    logging.info("Done.")
 
# ============================================================================
# 11. Entry point
# ============================================================================
 
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="EasyCORA-Landsat matchup processor")
    ap.add_argument('--verify', action='store_true',
                    help='Run verification pass (submitted -> exported/empty)')
    ap.add_argument('--csv-dir', type=Path, default=None,
                    help='Local folder containing downloaded CSVs (for --verify)')
    return ap.parse_args()
 
def main() -> None:
    setup_logging()
    args = parse_args()
 
    try:
        ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
    except Exception as e:
        logging.error(f"Earth Engine auth failed: {e}")
        sys.exit(1)
 
    if args.verify:
        verify(args.csv_dir)
    else:
        main_process()
 
if __name__ == '__main__':
    main()