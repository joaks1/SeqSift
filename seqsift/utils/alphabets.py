#! /usr/bin/env python

import copy

from seqsift.utils import iteritems

class ResidueAlphabet(object):
    def __init__(self,
            states,
            gap = None,
            ambiguity_codes = {}):

        self._states = None
        self._gap = None
        self._missing = None
        self._ambiguity_codes = {}
        self._reverse_ambiguity_codes = {}
        self._residue_ambiguity_codes = {}

        self.states = states
        self.gap = gap
        self.ambiguity_codes = ambiguity_codes

    def _set_states(self, states):
        self._states = self.standardize_states(states)

    def _get_states(self):
        return self._states

    states = property(_get_states, _set_states)

    def has_state(self, state):
        return (state in self.states)

    def _get_residues(self):
        return self.standardize_states(
                [s for s in self.states if s != self.gap])

    residues = property(_get_residues)

    def _set_gap_state(self, state):
        if state is None:
            self._gap = state
            return
        if not state in self.states:
            raise ValueError('Gap symbol {0!r} is not in state alphabet: '
                    '{1}'.format(state, self.states))
        if state in self.ambiguity_codes:
            raise ValueError('Gap state symbol {0!r} is already used as '
                    'ambiguity')
        self._gap = state

    def _get_gap_state(self):
        return self._gap

    gap = property(_get_gap_state, _set_gap_state)

    def _set_ambiguity_codes(self, symbol_to_states_dict):
        self._ambiguity_codes = {}
        self._reverse_ambiguity_codes = {}
        for ambig, states in iteritems(symbol_to_states_dict):
            for s in states:
                if not s in self.states:
                    raise ValueError('Ambiguity {0!r} maps to an invalid '
                            'state {0!r}'.format(ambig, s))
            tup = self.standardize_states(states)
            self._ambiguity_codes[ambig] = tup
            self._reverse_ambiguity_codes[tup] = ambig
            if tup == self.states:
                self._missing = ambig
            if not self.gap or (self.gap not in tup):
                self._residue_ambiguity_codes[ambig] = tup

    def _get_ambiguity_codes(self):
        return self._ambiguity_codes

    ambiguity_codes = property(_get_ambiguity_codes, _set_ambiguity_codes)

    def _get_missing_symbol(self):
        return self._missing

    missing = property(_get_missing_symbol)

    def _get_reverse_ambiguity_codes(self):
        return self._reverse_ambiguity_codes

    reverse_ambiguity_codes = property(_get_reverse_ambiguity_codes)

    def _get_residue_ambiguity_codes(self):
        return self._residue_ambiguity_codes

    residue_ambiguity_codes = property(_get_residue_ambiguity_codes)

    def _get_all_residue_codes(self):
        d = copy.deepcopy(self._residue_ambiguity_codes)
        for r in self.residues:
            d[r] = self.standardize_states([r])
        return d

    all_residue_codes = property(_get_all_residue_codes)

    def get_symbol(self, states):
        return self._reverse_ambiguity_codes.get(
                self.standardize_states(states),
                None)

    def standardize_states(self, states):
        return tuple(sorted(states))

    def get_valid_symbols(self):
        return list(set(list(self.states) + list(self.ambiguity_codes.keys())))

class DnaAlphabet(ResidueAlphabet):
    def __init__(self):
        ResidueAlphabet.__init__(self,
                states = ('A', 'C', 'G', 'T', '-'),
                gap = '-',
                ambiguity_codes = {
                        'R': ('A', 'G'),
                        'Y': ('C', 'T'),
                        'K': ('G', 'T'),
                        'M': ('A', 'C'),
                        'S': ('C', 'G'),
                        'W': ('A', 'T'),
                        'V': ('A', 'C', 'G'),
                        'H': ('A', 'C', 'T'),
                        'D': ('A', 'G', 'T'),
                        'B': ('C', 'G', 'T'),
                        'N': ('A', 'C', 'G', 'T'),
                        '?': ('A', 'C', 'G', 'T', '-'),
                        })

class RnaAlphabet(ResidueAlphabet):
    def __init__(self):
        ResidueAlphabet.__init__(self,
                states = ('A', 'C', 'G', 'U', '-'),
                gap = '-',
                ambiguity_codes = {
                        'R': ('A', 'G'),
                        'Y': ('C', 'U'),
                        'K': ('G', 'U'),
                        'M': ('A', 'C'),
                        'S': ('C', 'G'),
                        'W': ('A', 'U'),
                        'V': ('A', 'C', 'G'),
                        'H': ('A', 'C', 'U'),
                        'D': ('A', 'G', 'U'),
                        'B': ('C', 'G', 'U'),
                        'N': ('A', 'C', 'G', 'U'),
                        '?': ('A', 'C', 'G', 'U', '-'),
                        })

class ProteinAlphabet(ResidueAlphabet):
    def __init__(self):
        ResidueAlphabet.__init__(self,
                states = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*',
                    '-'),
                gap = '-',
                ambiguity_codes = {
                        'B': ('D', 'N'),
                        'Z': ('E', 'Q'),
                        'X': ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
                            'Y', '*'),
                        '?': ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
                            'Y', '*', '-'),
                        })

