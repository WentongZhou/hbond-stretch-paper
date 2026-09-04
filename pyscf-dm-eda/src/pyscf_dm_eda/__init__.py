"""Public API for the standalone PySCF density-matrix EDA package."""

from .eda import (
    Atom,
    EDAError,
    EDAResult,
    EDAValidationError,
    FragmentSpec,
    GridRunSummary,
    HARTREE_TO_KCAL_MOL,
    HARTREE_TO_EV,
    HARTREE_TO_KJ_MOL,
    IncompatibleFragmentError,
    PySCFEDA,
    PySCFGridRunner,
    SCFConfig,
    SCFConvergenceError,
    grid_main,
    read_xyz,
)

__version__ = "0.1.0"

__all__ = [
    "Atom",
    "EDAError",
    "EDAResult",
    "EDAValidationError",
    "FragmentSpec",
    "GridRunSummary",
    "HARTREE_TO_KCAL_MOL",
    "HARTREE_TO_EV",
    "HARTREE_TO_KJ_MOL",
    "IncompatibleFragmentError",
    "PySCFEDA",
    "PySCFGridRunner",
    "SCFConfig",
    "SCFConvergenceError",
    "grid_main",
    "read_xyz",
]
