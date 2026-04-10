"""
cht_cyclones package initialisation.

Exposes the primary public API: TropicalCyclone, TropicalCycloneTrack,
CycloneTrackDatabase, CycloneTrackDataset, and the track_selector helper.
"""

__version__ = "1.0.3"

# Import jtwc last — it depends on TropicalCyclone which must be defined first
import cht_cyclones.jtwc.jtwc as jtwc  # noqa: F401
from cht_cyclones.track import TropicalCycloneTrack
from cht_cyclones.track_database import CycloneTrackDatabase
from cht_cyclones.track_dataset import CycloneTrackDataset
from cht_cyclones.track_selector import track_selector
from cht_cyclones.tropical_cyclone_refactored import TropicalCyclone

__all__ = [
    "TropicalCycloneTrack",
    "CycloneTrackDatabase",
    "CycloneTrackDataset",
    "track_selector",
    "TropicalCyclone",
]
