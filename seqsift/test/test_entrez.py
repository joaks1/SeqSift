#! /usr/bin/env python

import os
import sys
import unittest
import itertools

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from seqsift.utils.entrez import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class testParseGiNumbersTestCase(unittest.TestCase):
    def test_empty_string(self):
        gis = parse_gi_numbers('')
        self.assertEqual(gis, [])
    
    def test_invalid_gis(self):
        gis = parse_gi_numbers('AA123, 221!#3, 123A34-123C344')
        self.assertEqual(gis, [])

    def test_invalid_range(self):
        gis = parse_gi_numbers('1234-1223')
        self.assertEqual(gis, [])
    
    def test_basic(self):
        gis = parse_gi_numbers('1212, 1234, 2300-2312')
        self.assertEqual(sorted(gis),
                ['1212', '1234'] + [str(x) for x in range(2300,2313)])

class testParseAccessionNumbers(unittest.TestCase):
    def test_empty_string(self):
        accs = parse_accession_numbers('')
        self.assertEqual(accs, [])

    def test_invalid_accessions(self):
        accs = parse_accession_numbers(
                '123434, AAAAAA344535, AU2342333A, AU1234')
        self.assertEqual(accs, [])
    
    def test_invalid_range(self):
        accs = parse_accession_numbers('AU23005-AU23004')
        self.assertEqual(accs, [])
        accs = parse_accession_numbers('AU23005-AX23006')
        self.assertEqual(accs, [])

    def test_basic(self):
        accs = parse_accession_numbers(
            'A12315, AAAAA256210, AAU23434-AAU23435, AU222222-222224')
        self.assertEqual(sorted(accs), sorted([
                'A12315', 'AAAAA256210', 'AAU23434', 'AAU23435',
                'AU222222', 'AU222223', 'AU222224']))

    
class GetEntrezDatabaseTestCase(unittest.TestCase):
    def test_error_handling(self):
        self.assertRaises(ValueError, get_entrez_database, 'bogus')

    def test_nuccore(self):
        db = get_entrez_database('dna')
        self.assertEqual(db, 'nuccore')
        db = get_entrez_database('rna')
        self.assertEqual(db, 'nuccore')

    def test_protein(self):
        db = get_entrez_database('aa')
        self.assertEqual(db, 'protein')
        db = get_entrez_database('protein')
        self.assertEqual(db, 'protein')

class ParseMixedGiListTestCase(unittest.TestCase):
    def setUp(self):
        self.id_list = ['354698776', '354698778', '354698780', '354698782']
        self.acc_list = ['JF314863', 'JF314864', 'JF314865', 'JF314866']

    def test_gi_numbers(self):
        gis, accs = parse_mixed_gi_list(self.id_list)
        self.assertEqual(accs, [])
        self.assertEqual(sorted(gis), sorted(self.id_list))

    def test_accession_numbers(self):
        gis, accs = parse_mixed_gi_list(self.acc_list)
        self.assertEqual(sorted(accs), sorted(self.acc_list))
        self.assertEqual(sorted(gis), [])

    def test_mixed_list(self):
        gis, accs = parse_mixed_gi_list(self.id_list[:2] + self.acc_list[2:])
        self.assertEqual(sorted(accs), sorted(self.acc_list[2:]))
        self.assertEqual(sorted(gis), sorted(self.id_list[:2]))

class GetPureGiNumbersTestCase(unittest.TestCase):
    def setUp(self):
        self.id_list = ['354698776', '354698778', '354698780', '354698782']
        self.acc_list = ['JF314863', 'JF314864', 'JF314865', 'JF314866']

    def test_gi_numbers(self):
        gis = get_pure_gi_numbers(self.id_list)
        self.assertEqual(sorted(gis), sorted(self.id_list))

    def test_accession_numbers(self):
        gis = get_pure_gi_numbers(self.acc_list)
        self.assertEqual(sorted(gis), sorted(self.id_list))

    def test_mixed_list(self):
        gis = get_pure_gi_numbers(self.id_list[:2] + self.acc_list[2:])
        self.assertEqual(sorted(gis), sorted(self.id_list))

    def test_mixed_list_overlapping(self):
        gis = get_pure_gi_numbers(self.id_list + self.acc_list)
        self.assertEqual(sorted(gis), sorted(self.id_list))

class GetGbHandleTestCase(SeqSiftTestCase):
    def setUp(self):
        self.id = '354698774'
        self.acc = 'JF314862'
        self.singleton_fasta = package_paths.data_path('JF314862.fasta')
        self.singleton_gb = package_paths.data_path('JF314862.gb')
        self.id_list = ['354698776', '354698778', '354698780', '354698782']
        self.acc_list = ['JF314863', 'JF314864', 'JF314865', 'JF314866']
        self.long_acc_list = [
                'JF314862',
                'JF314863',
                'JF314864',
                'JF314865',
                'JF314866',
                'JF314867',
                'JF314868',
                'JF314869',
                'JF314870',
                'JF314871',
                'JF314872',
                'JF314873',
                'JF314874',
                'JF314875',
                'JF314876',]
        self.ids = ','.join(self.id_list)
        self.multi_fasta = package_paths.data_path('JF314863-JF314866.fasta')
        self.multi_gb = package_paths.data_path('JF314863-JF314866.gb')
        self.long_multi_fasta = package_paths.data_path(
                'JF314862-JF314876.fasta')
        self.long_multi_gb = package_paths.data_path('JF314862-JF314876.gb')

    def test_singleton_fasta(self):
        h = get_gb_handle(self.id, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.singleton_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_singleton_gb(self):
        h = get_gb_handle(self.id, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.singleton_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()
    
    def test_singleton_accession(self):
        h = get_gb_handle(self.acc, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.singleton_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_fasta(self):
        h = get_gb_handle(self.ids, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_gb(self):
        h = get_gb_handle(self.ids, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_accession(self):
        h = get_gb_handle(self.acc_list, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_long_accession_gb(self):
        h = get_gb_handle(self.long_acc_list, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.long_multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_long_accession_fasta(self):
        h = get_gb_handle(self.long_acc_list, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.long_multi_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_singleton_fasta_tmp_file(self):
        h = get_gb_handle(self.id, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=True)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.singleton_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_singleton_gb_tmp_file(self):
        h = get_gb_handle(self.id, db='nuccore', rettype='gb',
            retmode='text', tmp_file=True)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.singleton_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_fasta_tmp_file(self):
        h = get_gb_handle(self.ids, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=True)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_gb_tmp_file(self):
        h = get_gb_handle(self.ids, db='nuccore', rettype='gb',
            retmode='text', tmp_file=True)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_fasta_id_list(self):
        h = get_gb_handle(self.id_list, db='nuccore', rettype='fasta',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='fasta', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_fasta, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

    def test_multi_gb_id_list(self):
        h = get_gb_handle(self.id_list, db='nuccore', rettype='gb',
            retmode='text', tmp_file=False)
        seqs1 = SeqIO.parse(h, format='gb', alphabet=IUPAC.ambiguous_dna)
        seqs2 = SeqIO.parse(self.multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)
        h.close()

class FetchGbSeqsTestCase(SeqSiftTestCase):
    def setUp(self):
        self.id_list = ['354698776', '354698778', '354698780', '354698782']
        self.multi_gb = package_paths.data_path('JF314863-JF314866.gb')

    def test_basic(self):
        seqs1 = fetch_gb_seqs(self.id_list, data_type='dna')
        seqs2 = SeqIO.parse(self.multi_gb, format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameData(seqs1, seqs2)

if __name__ == '__main__':
    unittest.main()
