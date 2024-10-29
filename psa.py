"""Global and local pairwise alignments between nucleotide/protein sequences.

The module uses needle/water from the EMBOSS package to compute an optimal
global/local alignment between a pair of sequences (query and subject).

Copyright 2022 Andrzej Zielezinski (a.zielezinski@gmail.com)
https://github.com/aziele/pairwise-sequence-alignment

Adapted to also use stretcher based on needle.
"""

from __future__ import annotations
from collections import namedtuple
from typing import Optional
import subprocess
import shutil
import random

__version__ = '1.0.1'

# Check whether needle is on PATH and marked as executable.
assert shutil.which('needle'), "needle not found (is emboss installed?)"


class PairwiseAlignment():
    """Object representing a pairwise ailgnment.

    Attributes:
        qid          Query sequence identifier
        sid          Subject sequence identifier
        qseq         Query unaligned sequence
        sseq         Subject unaligned sequence
        qaln         Query aligned sequence
        saln         Subject aligned sequence
        qstart       Start of alignment in query
        qend         End of alignment in query
        sstart       Start of alignment in subject
        send         End of alignment in subject
        length       Alignment length
        score        Alignment score
        nidentity    Number of identical matches in the alignment
        pidentity    Percentage of identical matches in the alignment
        nsimilarity  Number of positive-scoring matches in the alignment
        psimilarity  Percentage of positive-scoring matches in the alignment
        ngaps        Total number of gaps in the alignment
        pgaps        Total percentage of gaps in the alignment
        moltype      nucl/prot
        program      needle/water
        gapopen      Gap open penalty
        gapextend    Gap extension penalty
        matrix       Name of scoring matrix
        raw          Raw output obtained from EMBOSS' needle/water
    """

    def __init__(
        self,
        qid: str,
        sid: str,
        qseq: str, 
        sseq: str,
        qaln: str,
        saln: str,
        qstart: int,
        qend: int,
        sstart: int,
        send: int,
        length: int,
        score: float,
        nidentity: int,
        nsimilarity: int,
        ngaps: int,
        moltype: str,
        program: str, 
        gapopen: int,
        gapextend: int,
        matrix: str,
        raw: List[str]
    ):
        self.qid = qid
        self.sid = sid
        self.qseq = qseq
        self.sseq = sseq
        self.qaln = qaln
        self.saln = saln
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.length = length
        self.score = score
        self.nidentity = nidentity
        self.nsimilarity = nsimilarity
        self.ngaps = ngaps
        self.moltype = moltype
        self.program = program
        self.gapopen = gapopen
        self.gapextend = gapextend
        self.matrix = matrix
        self._raw = raw
        self.qlen = len(qseq)
        self.slen = len(sseq)

    @property
    def raw(self) -> str:
        """Raw output from EMBOSS' needle/water"""
        return "\n".join(self._raw)

    @property    
    def pidentity(self) -> float:
        """Percentage of identical matches"""
        return self.nidentity / self.length * 100

    @property    
    def psimilarity(self) -> float:
        """Percentage of positive-scoring matches"""
        return self.nsimilarity / self.length * 100

    @property    
    def pgaps(self) -> float:
        """Percentage of gaps"""
        return self.ngaps / self.length * 100

    def query_coverage(self) -> float:
        """Query coverage"""
        return (self.qend - self.qstart + 1) / self.qlen * 100

    def subject_coverage(self) -> float:
        """Query coverage"""
        return (self.send - self.sstart + 1) / self.slen * 100

    def fasta(self, wrap=70) -> str:
        """Returns pairwise alignment in FASTA/Pearson format."""
        lst = [
            f'>{self.qid} {self.qstart}-{self.qend}',
            *(self.qaln[i:i + wrap] for i in range(0, self.length, wrap)),
            f'>{self.sid} {self.sstart}-{self.send}',
            *(self.saln[i:i + wrap] for i in range(0, self.length, wrap)),
        ]
        return "\n".join(lst)

    def pvalue(self, n:int = 100) -> float:
        """Returns p-value of the alignment.

        The method shuffles a subject sequence many times (n=100) and calculates
        the alignment score between the query sequence and each shuffled subject
        sequence. It then counts how many times the alignment score was greater
        than or equal to the alignment score of original (unshuffled) query and
        subject sequences.

        Args:
            n: Number of times query is aligned to the shuffled subject sequence
        """
        slst = list(self.sseq)
        k = 0
        for i in range(n):
            random.shuffle(slst)
            sseq = "".join(slst)
            aln = align(
                program=self.program,
                moltype=self.moltype,
                qseq=self.qseq,
                sseq=sseq,
                qid='query',
                sid='subject',
                gapopen=self.gapopen,
                gapextend=self.gapextend,
                matrix=self.matrix
            )
            if aln.score >= self.score:
                k += 1
        return k / n

    def __str__(self) -> str:
        """Returns the raw EMBOSS output."""
        return self.raw

    def __len__(self) -> int:
        """Returns the alignment length."""
        return self.length

    def __iter__(self):
        """Iterates over the characters in the alignment."""
        return iter(zip(self.qaln, self.saln))

    def __getitem__(self, index):
        """Returns aligned characters at given position (index)."""
        return (self.qaln[index], self.saln[index])


SimpleAlignment = namedtuple(
    'SimpleAlignment', 
    ['qaln', 'saln', 'qstart', 'qend', 'sstart', 'send', 'length', 'score',
     'nidentity', 'nsimilarity', 'ngaps', 'output']
)

def align(
        program: Literal['needle', 'water'],
        moltype: Literal['prot', 'nucl'],
        qseq: str,
        sseq: str,
        qid: str = 'query',
        sid: str = 'subject',
        gapopen: int = 10,
        gapextend: int = 0.5,
        matrix: Optional[str] = None
        ) -> PairwiseAlignment:
    """Aligns two sequences, parses the output, and returns an alignment object.

    Args:
        program   An EMBOSS tool to run.
        moltype   A molecule type of query and subject sequences
        qseq      Query sequence
        sseq      Subject sequence
        qid       Query sequence identifier
        sid       Subject sequence identifier
        gapopen   Gap open penalty
        gapextend Gap extension penalty
        matrix    Name of scoring matrix

    Returns:
        A list of lines from EMBOSS output
    """
    if not matrix:
        matrix = 'EBLOSUM62' if moltype == 'prot' else 'EDNAFULL'
    handle = emboss_run(
        program=program,
        moltype=moltype,
        qseq=qseq,
        sseq=sseq,
        qid=qid,
        sid=sid,
        gapopen=gapopen,
        gapextend=gapextend,
        matrix=matrix
    )
    aln = emboss_parse(handle)
    return PairwiseAlignment(
        qid=qid,
        sid=sid,
        qaln=aln.qaln,
        saln=aln.saln,
        qseq=qseq,
        sseq=sseq,
        qstart=aln.qstart,
        qend=aln.qend,
        sstart=aln.sstart,
        send=aln.send,
        length=aln.length,
        score=aln.score,
        nidentity=aln.nidentity, 
        nsimilarity=aln.nsimilarity,
        ngaps=aln.ngaps,
        moltype=moltype,
        program=program,
        gapopen=gapopen,
        gapextend=gapextend,
        matrix=matrix,
        raw=aln.output
    )

def emboss_run(
        program: Literal['needle', 'water', 'stretcher'],
        moltype: Literal['prot', 'nucl'],
        qseq: str,
        sseq: str,
        qid: str,
        sid: str,
        gapopen: int,
        gapextend: int,
        matrix: str,
        ) -> Iterable[str]:
    """Aligns two sequences using EMBOSS and returns its output.

    Args:
        program   An EMBOSS tool to run.
        moltype   A molecule type of query and subject sequences
        qseq      Query sequence
        sseq      Subject sequence
        qid       Query sequence identifier
        sid       Subject sequence identifier
        gapopen   Gap open penalty
        gapextend Gap extension penalty
        matrix    Name of scoring matrix

    Returns:
        A list of lines from EMBOSS output

    Raises:
        CalledProcessError: if returncode of needle/water is non-zero
    """
    qualifiers = '-snucleotide1 -snucleotide2'
    if moltype == 'prot':
        qualifiers = '-sprotein1 -sprotein2'
    cmd = [
        f"{program} -stdout -auto",
        f"{qualifiers}",
        f"-asequence <(echo '>{qid}'; echo '{qseq}';)",
        f"-bsequence <(echo '>{sid}'; echo '{sseq}';)",
        f"-datafile {matrix}",
        f"-gapopen {gapopen}",
        f"-gapextend {gapextend}",
    ]
    process = subprocess.run(
        " ".join(cmd),
        shell=True,
        text=True,
        capture_output=True,
        executable='/bin/bash'
    )
    process.check_returncode()
    return process.stdout.splitlines()


def emboss_parse(handle: Iterable[str]) -> collections.namedtuple:
    """Parses EMBOSS output.

    Args:
        handle: A list of lines in EMBOSS output

    Returns:
        A namedtuple containing alignment data.
    """
    qaln = []
    saln = []
    output = []
    qpos = []
    spos = []
    is_output = False
    for line in handle:
        if line.startswith('#====='):
            is_output = True
        if is_output:
            output.append(line)

        # Skip blank/empty lines
        if not line.rstrip():
            continue

        is_header = True if line.startswith('#') else False
        cols = line.split()
        # Parse header
        if is_header:
            if line.startswith('# 1: '):
                qid = cols[2]
            elif line.startswith('# 2: '):
                sid = cols[2]
            elif line.startswith('# Length:'):
                length = int(cols[2])
            elif line.startswith('# Identity: '):
                nidentity = int(cols[2].split('/')[0])
            elif line.startswith('# Similarity: '):
                nsimilarity = int(cols[2].split('/')[0])
            elif line.startswith('# Gaps: '):
                ngaps = int(cols[2].split('/')[0])
            elif line.startswith('# Score: '):
                score = float(cols[2])
        # Parse alignment
        else:
            if line.startswith(qid):
                qaln.append(cols[2])
                qpos.append(int(cols[1]))
                qpos.append(int(cols[3]))
            elif line.startswith(sid):
                saln.append(cols[2])
                spos.append(int(cols[1]))
                spos.append(int(cols[3]))

    return SimpleAlignment(
        "".join(qaln),
        "".join(saln),
        min(qpos),
        max(qpos),
        min(spos),
        max(spos),
        length,
        score,
        nidentity,
        nsimilarity,
        ngaps,
        output
    )


def water(*args, **kwargs):
    """Aligns two sequences using EMBOSS water."""
    return align('water', *args, **kwargs)


def needle(*args, **kwargs):
    """Aligns two sequences using EMBOSS needle."""
    return align('needle', *args, **kwargs)


def stretcher(*args, **kwargs):
    """Aligns two sequences using EMBOSS stretcher."""
    return align('stretcher', *args, **kwargs)
