"""Models for storing BAM/CRAM statistics information for targeted or whole genome sequencing data.

This includes statistics of the coverage as well as basic read alignment statistics such as
duplicate rate.
"""

import enum


class SequencingType(enum.Enum):
    """The sequncing types."""

    #: Whole-genome sequencing.
    WGS = "wgs"
    #: Targeted sequencing (exome/panel).
    TARGETED = "targeted"
    # TODO: RNA-seq


class BasicAlignmentStats:
    """Basic alignment statistics."""


class TargetedCoverageStats:
    """Coverage statistics for targeted sequencing data."""


class WgsCoverageStats:
    """Coverage statistics for WGS data."""
