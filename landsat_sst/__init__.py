"""
python_landsat_sst — Earth Engine Python API module for Arctic Landsat SST retrieval.

Implements the interaction-term algorithm derived in
python_scripts_to_derive_Algorithm/derive_SST-COEFFS-main.py.

Quick start::

    import ee
    from python_landsat_sst import get_collection
    from python_landsat_sst.water_mask import build_water_mask

    ee.Initialize(project='your-project')

    roi = ee.Geometry.Rectangle([-53, 70.3, -49.5, 71.1])
    col = get_collection('L8', '2022-07-01', '2022-08-01', roi)

    image = col.first()
    qa    = image.select('QA_PIXEL')
    sst   = image.select('SST').updateMask(build_water_mask(qa, 'L8'))
"""

from .collection import get_collection, COLLECTION_INFO, RGB_BANDS
from .water_mask  import build_water_mask, apply_water_mask
from .algorithm   import add_sst_band
from .params      import COEFFS, SCALER_PARAMS, TIR_BANDS, SENSORS

__all__ = [
    'get_collection',
    'build_water_mask',
    'apply_water_mask',
    'add_sst_band',
    'COEFFS',
    'SCALER_PARAMS',
    'TIR_BANDS',
    'SENSORS',
    'COLLECTION_INFO',
    'RGB_BANDS',
]
