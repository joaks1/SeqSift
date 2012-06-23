import sys
import os
import copy

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.utils.dataio import get_seq_iter
from seqsift.utils.entrez import get_gb_handle
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class Fragment(SeqRecord):
    def __init__(self, 
                 seq_record,
                 start_site,
                 end_site,
                 overhang,
                 five_prime_terminus=False,
                 three_prime_terminus=False,
                 **kwargs):
        if start_site > end_site:
            raise InvalidFragmentError(
                "A Fragment's start position ({0}) cannot be greater than "
                "its end position ({1})".format(start_site, end_site))
        if isinstance(seq_record, SeqRecord):
            SeqRecord.__init__(self,
                    seq = seq_record.seq,
                    id = seq_record.id,
                    name = seq_record.name,
                    description = seq_record.description,
                    letter_annotations = seq_record.letter_annotations,
                    annotations = seq_record.annotations,
                    features = seq_record.features,
                    dbxrefs = seq_record.dbxrefs)
        else:
            SeqRecord.__init__(self,
                    seq = Seq(str(seq_record).upper(), 
                            alphabet=IUPAC.ambiguous_dna),
                    id = kwargs.get('id', '<unknown id>'),
                    name = kwargs.get('name', '<unknown name>'),
                    description = kwargs.get('description',
                            '<unknown description>'),
                    letter_annotations = kwargs.get('letter_annotations',
                            None),
                    annotations = kwargs.get('annotations', None),
                    features = kwargs.get('features', None),
                    dbxrefs = kwargs.get('dbxrefs', None))
        self.start = int(start_site)
        self.end = int(end_site)
        self.length = self.end - self.start + 1
        if self.length != len(self):
            raise InvalidFragmentError(
                "The fragment's sequence length ({0}) does not match its "
                "start ({1}) and end ({2}) sites".format(len(self),
                        self.start, self.end))
        self.overhang = int(overhang)
        self.five_prime_terminus = five_prime_terminus
        self.three_prime_terminus = three_prime_terminus

class InvalidFragmentError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class RecognitionSeq(SeqRecord):
    """
    A recognition sequence of a restriction enzyme. The class currently only
    supports unambiguous DNA recognition sequences.
    """
    def __init__(self, seq_record, cut_site, **kwargs):
        """
        __init__ requires a string representation of the recognition
        sequence, `seq_str`, and the one-based index of the last base
        before the cut site `cut_site`.

        E.g., NotI:
            5'---GC     GGCCGC---3'
            3'---CGCCGG     CG---5'
            noti = RecognitionSeq('GCGGCCGC', 2)

        E.g., EcoRV:
            5'---GAT  ATC---3'
            3'---CTA  TAG---5'
            ecorv = RecognitionSeq('GATATC', 3)
        """
        if isinstance(seq_record, SeqRecord):
            SeqRecord.__init__(self,
                    seq = seq_record.seq,
                    id = seq_record.id,
                    name = seq_record.name,
                    description = seq_record.description,
                    letter_annotations = seq_record.letter_annotations,
                    annotations = seq_record.annotations,
                    features = seq_record.features,
                    dbxrefs = seq_record.dbxrefs)
        else:
            SeqRecord.__init__(self,
                    seq = Seq(str(seq_record).upper(), 
                            alphabet=IUPAC.unambiguous_dna),
                    id = kwargs.get('id', '<unknown id>'),
                    name = kwargs.get('name', '<unknown name>'),
                    description = kwargs.get('description',
                            '<unknown description>'),
                    letter_annotations = kwargs.get('letter_annotations',
                            None),
                    annotations = kwargs.get('annotations', None),
                    features = kwargs.get('features', None),
                    dbxrefs = kwargs.get('dbxrefs', None))
        for base in str(self.seq):
            if base.upper() not in ['A', 'C', 'G', 'T']:
                raise InvalidRecognitionSeqError(
                    "Invalid base {0!r}:\n\t".format(base.upper()) + 
                    "RecognitionSeq only supports unambiguous DNA")
        self.cut_site = int(cut_site)
        if self.cut_site < 0 or self.cut_site > len(self):
            raise InvalidRecognitionSeqError(
                "Invalid cut_site {0} for recognition sequence {1}".format(
                        self.cut_site, str(self.seq)))
        self.overhang = len(self) - (2*self.cut_site)

    def _potential_overlap(self):
        for i in range(1, len(self)):
            if ((str(self.seq)[:i] == str(self.seq)[-i:]) and
                    ((self.cut_site < i) or
                     (self.cut_site > len(self) - i))):
                return True
        return False

    potential_overlap = property(_potential_overlap)

    def _palindrome(self):
        if str(self.seq) == str(self.seq.reverse_complement()):
            return True
        return False

    palindrome = property(_palindrome)

    def digest(self, seq_record):
        if not self.palindrome:
            _LOG.warning("WARNING: "
                "The recognition sequence {0!r} is not a palindrome.\n"
                "The digest method will only perform a single-stranded "
                "digest of the DNA molecule.\nThe result might not be "
                "very meaningful.".format(str(self.seq)))
        if self.potential_overlap:
            _LOG.warning("WARNING: "
                "The recognition sequence {0!r} has the potential for "
                "overlapping recognition sites.\nThe resulting digest "
                "will only represent one of the possibly many "
                "potential outcomes.".format(str(self.seq)))
        start_site = 1
        cut_buffer = 0
        five_prime_terminus = True
        three_prime_terminus = False
        for i in range(0, len(seq_record.seq) - len(self) + 1):
            if cut_buffer > 0:
                cut_buffer -= 1
                continue
            if (str(seq_record.seq[i: i + len(self)]).upper() == str(self.seq)):
                if i == 0 and self.cut_site == 0:
                    continue
                if start_site != 1:
                    five_prime_terminus = False
                if (i + self.cut_site + 1) > len(seq_record.seq):
                    three_prime_terminus = True
                    if five_prime_terminus or self.overhang >= 0:
                        ohang = 0
                    else:
                        ohang = abs(self.overhang)
                    yield Fragment(
                            seq_record = seq_record[start_site - 1:],
                            start_site = start_site,
                            end_site = len(seq_record.seq),
                            overhang = ohang,
                            five_prime_terminus = five_prime_terminus,
                            three_prime_terminus = three_prime_terminus)
                    break
                end_site = i + self.cut_site
                if five_prime_terminus and self.overhang < 0:
                    ohang = 0
                else:
                    ohang = abs(self.overhang)
                yield Fragment(
                        seq_record = seq_record[start_site - 1: end_site],
                        start_site = start_site,
                        end_site = end_site,
                        overhang = ohang,
                        five_prime_terminus = five_prime_terminus,
                        three_prime_terminus = three_prime_terminus)
                start_site = i + self.cut_site + 1
                cut_buffer = self.cut_site - 1
                five_prime_terminus = False
        if not three_prime_terminus:
            if five_prime_terminus or self.overhang >= 0:
                ohang = 0
            else:
                ohang = abs(self.overhang)
            yield Fragment(seq_record = seq_record[start_site - 1:],
                           start_site = start_site,
                           end_site = len(seq_record.seq),
                           overhang = ohang,
                           five_prime_terminus = five_prime_terminus,
                           three_prime_terminus = True)

class InvalidRecognitionSeqError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class DigestSummary(object):
    def __init__(self, recognition_seq, seq_record, extra_length=0,
            include_overhang=True):
        if isinstance(recognition_seq, str):
            self.recognition_seq = recognition_seq.upper()
            rs = RecognitionSeq(recognition_seq)
        else:
            rs = recognition_seq
            self.recognition_seq = str(recognition_seq.seq).upper()
        self.molecule_id = seq_record.id
        self.molecule_name = seq_record.name
        self.molecule_description = seq_record.description
        self.molecule_length = len(seq_record.seq)
        self.length_distribution = {}
        for fragment in rs.digest(seq_record):
            if include_overhang:
                l = fragment.length + fragment.overhang + extra_length
            else:
                l = fragment.length + extra_length
            if l not in self.length_distribution.keys():
                self.length_distribution[l] = 0
            self.length_distribution[l] += 1
