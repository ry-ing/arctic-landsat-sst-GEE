"""
Main Landsat SST collection builder 

"""

import ee
from .algorithm import add_sst_band
from .merra import add_tcwv_band

# Collection-2 Tier-1 
COLLECTION_INFO = {
    'L4': {
        'TOA': 'LANDSAT/LT04/C02/T1_TOA',
        'SR':  'LANDSAT/LT04/C02/T1_L2',
        'TIR': ['B6'],
        'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
    },
    'L5': {
        'TOA': 'LANDSAT/LT05/C02/T1_TOA',
        'SR':  'LANDSAT/LT05/C02/T1_L2',
        'TIR': ['B6'],
        'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
    },
    'L7': {
        'TOA': 'LANDSAT/LE07/C02/T1_TOA',
        'SR':  'LANDSAT/LE07/C02/T1_L2',
        'TIR': ['B6_VCID_1'],   # B6_VCID_2 excluded; QA_PIXEL comes from SR
        'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
    },
    'L8': {
        'TOA': 'LANDSAT/LC08/C02/T1_TOA',
        'SR':  'LANDSAT/LC08/C02/T1_L2',
        'TIR': ['B10'],         
        'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
    },
    'L9': {
        'TOA': 'LANDSAT/LC09/C02/T1_TOA',
        'SR':  'LANDSAT/LC09/C02/T1_L2',
        'TIR': ['B10'],
        'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
    },
}


RGB_BANDS = {
    'L4': ['SR_B3', 'SR_B2', 'SR_B1'],
    'L5': ['SR_B3', 'SR_B2', 'SR_B1'],
    'L7': ['SR_B3', 'SR_B2', 'SR_B1'],
    'L8': ['SR_B4', 'SR_B3', 'SR_B2'],
    'L9': ['SR_B4', 'SR_B3', 'SR_B2'],
}


def _scale_sr(image):
    """Apply Collection-2 Level-2 SR scale factors to optical bands only."""
    optical = image.select('SR_B.*')
    scaled = optical.multiply(0.0000275).add(-0.2)
    return image.addBands(scaled, None, True)   # overwrite originals in-place


def get_collection(sensor, start_date, end_date, geometry):
    """
    Build an SST ImageCollection for the given sensor and filter.

    The returned collection has these bands on every image:
      SST    — sea surface temperature (°C), all pixels
      bt_K   — brightness temperature (K) of the primary thermal band
      TCWV   — total column water vapour (kg m⁻²) from MERRA-2
      QA_PIXEL — raw Landsat quality bitmask (for masking downstream)
      SR_B*  — scaled surface reflectance (reflectance units, 0–1)

    Returns:
        ee.ImageCollection
    """
    if sensor not in COLLECTION_INFO:
        raise ValueError(f"sensor must be one of {list(COLLECTION_INFO)}; got '{sensor}'")

    info = COLLECTION_INFO[sensor]

    landsat_sr = (
        ee.ImageCollection(info['SR'])
          .filterDate(start_date, end_date)
          .filterBounds(geometry)
          .map(add_tcwv_band)
          .map(_scale_sr)
    )

    landsat_toa = (
        ee.ImageCollection(info['TOA'])
          .filterDate(start_date, end_date)
          .filterBounds(geometry)
    )

    vis_bands = info['SR_BANDS'] + ['TCWV']
    tir_bands = info['TIR']

    # Combine SR (optical + QA + TCWV) with TOA (thermal band) by system:index
    merged = landsat_sr.select(vis_bands).combine(landsat_toa.select(tir_bands), True)

    return merged.map(add_sst_band(sensor))
