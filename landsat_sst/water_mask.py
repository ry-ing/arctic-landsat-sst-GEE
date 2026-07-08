"""
QA_PIXEL-based water masking for Landsat Collection 2.
"""

import ee


def build_water_mask(qa, sensor):
    """
    Return a binary mask (1 = valid water pixel, 0 = masked) from QA_PIXEL.

    Keeps only pixels that are:
      - water (bit 7 = 1)
      - AND free of fill, dilated cloud, cloud, cloud shadow, snow
      - AND cloud/shadow/snow confidence < 'medium' (bits 8-13 < 2)
      - AND for L8/L9: free of cirrus (bit 2, bits 14-15)

    """
    is_water = qa.bitwiseAnd(1 << 7).neq(0)

    good_qa = (
        qa.bitwiseAnd(1 << 0).eq(0)               # no fill
          .And(qa.bitwiseAnd(1 << 1).eq(0))        # no dilated cloud
          .And(qa.bitwiseAnd(1 << 3).eq(0))        # no cloud
          .And(qa.bitwiseAnd(1 << 4).eq(0))        # no cloud shadow
          .And(qa.bitwiseAnd(1 << 5).eq(0))        # no snow/ice
          .And(qa.rightShift(8).bitwiseAnd(3).lt(2))    # cloud confidence < 2
          .And(qa.rightShift(10).bitwiseAnd(3).lt(2))   # shadow confidence < 2
          .And(qa.rightShift(12).bitwiseAnd(3).lt(2))   # snow confidence < 2
    )

    if sensor in ('L8', 'L9'):
        good_qa = (
            good_qa
            .And(qa.bitwiseAnd(1 << 2).eq(0))              # no cirrus
            .And(qa.rightShift(14).bitwiseAnd(3).lt(2))    # cirrus confidence < 2
        )

    return is_water.And(good_qa)


def apply_water_mask(image, sensor):
    qa = image.select('QA_PIXEL')
    return image.updateMask(build_water_mask(qa, sensor))
