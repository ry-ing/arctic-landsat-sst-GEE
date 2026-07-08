"""
Landsat SST retrieval 
Algorithm:
  Tb_scaled   = (Tb_K  - MEDIAN[0]) / IQR[0]    [RobustScaler on Tb]
  TCWV_scaled = (TCWV  - MEDIAN[1]) / IQR[1]    [RobustScaler on TCWV]
  SST_K       = A * Tb_scaled + B * (Tb_scaled * TCWV_scaled) + Intercept + D
  SST_C       = SST_K - 273.15
"""

import ee
from .params import COEFFS, SCALER_PARAMS, TIR_BANDS


def add_sst_band(sensor):
    """
    Return a map function that appends 'SST' (°C) and 'bt_K' (K) bands.
    Returns:
        Callable[[ee.Image], ee.Image]
    """
    coeffs = COEFFS[sensor]
    scaler = SCALER_PARAMS[sensor]
    median_tb,   median_tcwv = scaler['MEDIAN']
    iqr_tb,      iqr_tcwv   = scaler['IQR']
    tir_band = TIR_BANDS[sensor]

    a     = coeffs['A_Tb']
    b     = coeffs['B_Interaction']
    inter = coeffs['Intercept']
    d     = coeffs['D']

    def _mapper(image):
        tb   = image.select(tir_band).rename('Tb')
        tcwv = image.select('TCWV')

        # RobustScaler: z = (x - median) / IQR
        tb_scaled   = tb.subtract(median_tb).divide(iqr_tb)
        tcwv_scaled = tcwv.subtract(median_tcwv).divide(iqr_tcwv)

        interaction = tb_scaled.multiply(tcwv_scaled)
        sst_kelvin  = (tb_scaled.multiply(a)
                                .add(interaction.multiply(b))
                                .add(inter)
                                .add(d))

        sst_celsius = sst_kelvin.subtract(273.15).rename('SST')

        return image.addBands(sst_celsius).addBands(tb.rename('bt_K'))

    return _mapper
