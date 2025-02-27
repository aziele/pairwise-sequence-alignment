#!/usr/bin/env python3

import shutil
import unittest

import psa

class Test(unittest.TestCase):

    def test_is_emboss_installed(self):
        self.assertIsNotNone(shutil.which('needle'))
        self.assertIsNotNone(shutil.which('water'))
        self.assertIsNotNone(shutil.which('stretcher'))

    def test_needle_dna(self):
        aln = psa.needle(moltype='nucl', qseq='ATGCTAGATA', sseq='ATGCTAGTTA')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)

    def test_water_dna(self):
        aln = psa.water(moltype='nucl', qseq='ATGCTAGTTA', sseq='ATCCT')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)   

    def test_stretcher_dna(self):
        aln = psa.stretcher(moltype='nucl', qseq='ATGCTAGATA', sseq='ATGCTAGTTA')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)   
        
    def test_needle_protein(self):
        aln = psa.needle(moltype='prot', qseq='MERILIIMTGG', sseq='MEKILILM')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)

    def test_water_protein(self):
        aln = psa.water(moltype='prot', qseq='MERI', sseq='MEKILILM')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)

    def test_stretcher_protein(self):
        aln = psa.stretcher(moltype='prot', qseq='MERILIIMTGG', sseq='MEKILILM')
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)

    def test_needle_dna_qid_sid_1(self):
        qid = 'dna1'
        sid = 'dna2'
        qseq = 'ATGCTAGATA'
        sseq = 'ATGCTAGTTA'
        aln = psa.needle(moltype='nucl', qseq=qseq, sseq=sseq, qid=qid, sid=sid)
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)
        self.assertEqual(aln.qid, qid)
        self.assertEqual(aln.sid, sid)

    def test_needle_dna_qid_sid_2(self):
        qid = 'query1'
        sid = 'query10'
        qseq = 'ATGCTAGATA'
        sseq = 'ATGCTAGTTA'
        aln = psa.needle(moltype='nucl', qseq=qseq, sseq=sseq, qid=qid, sid=sid)
        self.assertEqual(len(aln.qaln), len(aln.saln))
        self.assertEqual(len(aln.qaln), aln.length)
        self.assertEqual(aln.qid, qid)
        self.assertEqual(aln.sid, sid)


class TestPairwiseAlignment(unittest.TestCase):

    def setUp(self):
        self.aln1 = psa.PairwiseAlignment(
            qid='query',
            sid='subject',
            qseq='MERILIIMTGGTITSIRDDEVTLELLIDYRKRFGDTQRFDIVKLMNIIS',
            sseq='MERILIIMTGGTISSIKKENILNVDDEVTLELLIDYRKRFGDSQKFDIVILIIS',
            qaln='MERILIIMTGGTITSIR-------DDEVTLELLIDYRKRFGDTQRFDIVKLMNIIS',
            saln='MERILIIMTGGTISSIKKENILNVDDEVTLELLIDYRKRFGDSQKFDIVIL--IIS',
            qstart=1, 
            qend=49,
            sstart=1,
            send=54,
            length=56,
            score=183.5,
            nidentity=42,
            nsimilarity=46,
            ngaps=9,
            moltype='prot',
            program='needle',
            gapopen=10,
            gapextend=0.5,
            matrix='EBLOSUM62',
            raw=[]  # Skipped for code readability.
        )

        self.aln2 = psa.PairwiseAlignment(
            qid='oddly_long_sequence_identifier1',
            sid='oddly_long_sequence_identifier2',
            qseq='ATGCTAGTAGTTGATTTTTT',
            sseq='ATGCTAGTAGATGAT',
            qaln='ATGCTAGTAGTTGAT',
            saln='ATGCTAGTAGATGAT',
            qstart=1,
            qend=15,
            sstart=1,
            send=15,
            length=15,
            score=66.0,
            nidentity=14,
            nsimilarity=14,
            ngaps=0,
            moltype='nucl',
            program='water',
            gapopen=10,
            gapextend=0.5,
            matrix='EDNAFULL',
            raw=[]  # Skipped for code readability.
        )

    def test_pidentity(self):
        self.assertEqual(self.aln1.pidentity, 75.0)
        self.assertAlmostEqual(self.aln2.pidentity, 93.3333333)

    def test_psimilarity(self):
        self.assertAlmostEqual(self.aln1.psimilarity, 82.1428571)
        self.assertAlmostEqual(self.aln2.psimilarity, 93.3333333)

    def test_pgaps(self):
        self.assertAlmostEqual(self.aln1.pgaps, 16.0714286)
        self.assertAlmostEqual(self.aln2.pgaps, 0)

    def test_query_coverage(self):
        self.assertAlmostEqual(self.aln1.query_coverage(), 100.0)
        self.assertAlmostEqual(self.aln2.query_coverage(), 75.0)

    def test_subject_coverage(self):
        self.assertAlmostEqual(self.aln1.subject_coverage(), 100.0)
        self.assertAlmostEqual(self.aln2.subject_coverage(), 100.0)

    def test_len(self):
        self.assertEqual(len(self.aln1), 56)
        self.assertEqual(len(self.aln2), 15)

    def test_getitem(self):
        self.assertEqual(self.aln1[13], ('T', 'S'))
        self.assertEqual(self.aln2[1], ('T', 'T'))

    def test_iter(self):
        self.assertEqual(list(iter(self.aln1))[0], ('M', 'M'))
        self.assertEqual(list(iter(self.aln2))[0], ('A', 'A'))

    def test_fasta_aln1(self):
        fasta = [
            '>query 1-49',
            'MERILIIMTGGTITSIR-------DDEVTLELLIDYRKRFGDTQRFDIVKLMNIIS',
            '>subject 1-54',
            'MERILIIMTGGTISSIKKENILNVDDEVTLELLIDYRKRFGDSQKFDIVIL--IIS',
        ]
        self.assertEqual(self.aln1.fasta(), "\n".join(fasta))

    def test_fasta_aln2(self):
        fasta = [
            '>oddly_long_sequence_identifier1 1-15',
            'ATGCTAGTAGTTGAT',
            '>oddly_long_sequence_identifier2 1-15',
            'ATGCTAGTAGATGAT',
        ]
        self.assertEqual(self.aln2.fasta(), "\n".join(fasta))

    def test_fasta_aln1_wrap10(self):
        fasta = [
            '>query 1-49',
            'MERILIIMTG',
            'GTITSIR---',
            '----DDEVTL',
            'ELLIDYRKRF',
            'GDTQRFDIVK',
            'LMNIIS',
            '>subject 1-54',
            'MERILIIMTG',
            'GTISSIKKEN',
            'ILNVDDEVTL',
            'ELLIDYRKRF',
            'GDSQKFDIVI',
            'L--IIS',
        ]
        self.assertEqual(self.aln1.fasta(wrap=10), "\n".join(fasta))

    def test_fasta_aln2_wrap10(self):
        fasta = [
            '>oddly_long_sequence_identifier1 1-15',
            'ATGCTAGTAG',
            'TTGAT',
            '>oddly_long_sequence_identifier2 1-15',
            'ATGCTAGTAG',
            'ATGAT',
        ]
        self.assertEqual(self.aln2.fasta(wrap=10), "\n".join(fasta))

unittest.main()
