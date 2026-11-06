"""Data models and validation for gfapdefinisher policy parameters."""

from dataclasses import dataclass
from typing import Generic, TypeVar

T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class PolicyValue(Generic[T]):
    """A parameter value guarded with required justification and provenance.

    Zero arbitrary magic numbers: every parameter in gfapdefinisher must document
    its justification and origin. If a value is uncalibrated, it must be explicitly
    flagged as provisional.
    """

    value: T
    justification: str
    provenance: str
    provisional: bool = False

    def __post_init__(self) -> None:
        if not self.justification or not self.justification.strip():
            raise ValueError(f"Policy parameter {self.value!r} must have a non-empty justification.")
        if not self.provenance or not self.provenance.strip():
            raise ValueError(f"Policy parameter {self.value!r} must have a non-empty provenance.")
