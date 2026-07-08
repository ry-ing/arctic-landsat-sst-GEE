"""
MERRA-2 total column water vapour (TCWV)

Adds the closest-in-time MERRA-2 TQV band to a Landsat image within a ±1-hour
window, bilinearly resampled to the Landsat 30 m grid.  Returns a constant
-999 image if no MERRA-2 granule is found 
"""

import ee

_MERRA_COLLECTION = 'NASA/GSFC/MERRA/slv/2'
_MERRA_BAND = 'TQV'


def add_tcwv_band(image):
    """
    Map function: adds a 'TCWV' band (kg m⁻²) to a Landsat image.

    Finds the MERRA-2 image closest in time to the Landsat overpass within
    ±1 hour, reprojects it to the Landsat scene, and adds it as
    the 'TCWV' band.

    Returns:
        ee.Image — input image with 'TCWV' band added.
    """
    landsat_time = image.date()

    def _set_time_diff(merra_image):
        diff = (ee.Number(merra_image.date().millis())
                  .subtract(landsat_time.millis())
                  .abs())
        return merra_image.set('DateDist', diff)

    merra_window = (
        ee.ImageCollection(_MERRA_COLLECTION)
          .filterDate(landsat_time.advance(-1, 'hour'),
                      landsat_time.advance(1, 'hour'))
          .select(_MERRA_BAND)
    )

    merra_sorted = merra_window.map(_set_time_diff).sort('DateDist')
    merra_size = merra_sorted.size()

    merra_closest = ee.Image(
        ee.Algorithms.If(
            merra_size.eq(0),
            ee.Image.constant(-999.0).rename(_MERRA_BAND),
            merra_sorted.first()
        )
    )

    target_proj = image.select(0).projection()
    tcwv = (
        merra_closest.select(_MERRA_BAND)
                     .resample('bilinear')
                     .reproject(crs=target_proj)
                     .rename('TCWV')
    )

    return image.addBands(tcwv, ['TCWV'])
