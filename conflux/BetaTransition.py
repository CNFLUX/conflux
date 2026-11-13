# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# SPDX-License-Identifier: MIT

"""Decode database transition assignments without overloading BSG's spin change."""

from operator import index
import re
from typing import NamedTuple, Optional


class Transition(NamedTuple):
    deltaJ: Optional[int]
    parity_change: Optional[bool]
    forbidden_order: Optional[int]
    unique: Optional[bool]

    @property
    def shape_model(self):
        if self.forbidden_order == 0:
            return 'allowed'
        if self.unique:
            return 'unique_forbidden'
        return 'allowed_approximation'

    @property
    def spectral_L(self):
        return self.deltaJ if self.shape_model == 'unique_forbidden' else 0


UNKNOWN_TRANSITION = Transition(None, None, None, None)


def decode_djpi(djpi):
    """Return the selected transition, all candidates, and an ambiguity flag.

    A trailing '-' denotes parity change. Select the least forbidden known
    assignment, preferring nonunique transitions at equal order. An unknown
    candidate leaves the assignment unresolved and uses an allowed approximation.
    """
    candidates = []
    for token in str(djpi).split(','):
        match = re.fullmatch(r'(\d+)(-?)', token.strip())
        if match is None:
            candidates.append(UNKNOWN_TRANSITION)
            continue
        deltaJ = int(match[1])
        parity_change = bool(match[2])
        order = max(deltaJ - 1, 0)
        if bool(order % 2) != parity_change:
            order += 1
        unique = order > 0 and deltaJ == order + 1
        candidates.append(Transition(deltaJ, parity_change, order, unique))

    distinct = set(candidates)
    if UNKNOWN_TRANSITION in distinct:
        selected = UNKNOWN_TRANSITION
    else:
        selected = min(distinct, key=lambda t: (t.forbidden_order, t.unique, t.deltaJ))
    return selected, tuple(candidates), len(distinct) > 1


class BranchTransitionMixin:
    """Keep transition metadata synchronized when branches are edited."""

    @property
    def dJpi(self):
        return self._dJpi

    @dJpi.setter
    def dJpi(self, assignment):
        if not isinstance(assignment, str):
            raise TypeError('dJpi must be a string such as "0", "1-", or "2-"')
        transition, candidates, ambiguous = decode_djpi(assignment)
        self._dJpi = assignment
        self._transition = transition
        self.transition_candidates = candidates
        self.transition_ambiguous = ambiguous
        if hasattr(self, 'Parameters'):
            self.UpdateParams()

    @property
    def deltaJ(self):
        return self._transition.deltaJ

    @property
    def parity_change(self):
        return self._transition.parity_change

    @property
    def forbidden_order(self):
        return self._transition.forbidden_order

    @property
    def unique(self):
        return self._transition.unique

    @property
    def shape_model(self):
        if self.custom_func is not None:
            return 'custom'
        return self._transition.shape_model

    def transition_parameters(self):
        return {
            'L': self._transition.spectral_L,
            'deltaJ': self.deltaJ,
            'parity_change': self.parity_change,
            'forbidden_order': self.forbidden_order,
            'unique': self.unique,
            'shape_model': self.shape_model,
            'dJpi': self.dJpi,
            'transition_ambiguous': self.transition_ambiguous,
        }


def validate_spectral_L(L):
    """BSG accepts zero for allowed shapes or deltaJ >= 2 for unique forbidden."""
    L = index(L)
    if L < 0 or L == 1:
        raise ValueError('Spectral L must be 0 (allowed) or deltaJ >= 2 (unique forbidden)')
    return L
