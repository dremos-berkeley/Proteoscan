"""Decoy database construction in the PEAKS "decoy fusion" spirit.

PEAKS DB (Zhang et al. 2012) appends a shuffled copy of each protein to the protein itself,
separated by 'J', so target and decoy halves share one entry and receive the same protein-level
bonuses. Search engines we drive (Sage, Comet) need separate entries, so we emit, for every
target protein, a decoy entry with the same accession prefixed by DECOY_TAG whose sequence is a
random shuffle of the target (the shuffle keeps the residue composition and length). FDR is then
computed as decoys / targets exactly as PEAKS does.
"""
from __future__ import annotations

import hashlib
import random
from pathlib import Path
from typing import Iterator

DECOY_TAG = "#DECOY#"


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def shuffle_protein(seq: str, seed: int, keep_cleavage_sites: bool = False, cleave_at: str = "KR") -> str:
    """Shuffle residues. With keep_cleavage_sites the K/R positions stay fixed so that the decoy
    tryptic peptide length distribution matches the target (PEAKS's paper describes a plain
    shuffle; keep_cleavage_sites=False reproduces that)."""
    rng = random.Random(seed)
    if not keep_cleavage_sites:
        chars = list(seq)
        rng.shuffle(chars)
        return "".join(chars)
    movable_idx = [i for i, c in enumerate(seq) if c not in cleave_at]
    movable = [seq[i] for i in movable_idx]
    rng.shuffle(movable)
    out = list(seq)
    for i, c in zip(movable_idx, movable):
        out[i] = c
    return "".join(out)


def build_target_decoy(fasta_in: str | Path, fasta_out: str | Path, seed: int = 0,
                       keep_cleavage_sites: bool = False, decoy_tag: str = DECOY_TAG) -> tuple[int, int]:
    """Write targets followed by shuffled decoys. Returns (n_targets, n_decoys)."""
    n = 0
    entries = list(read_fasta(fasta_in))
    with open(fasta_out, "w") as out:
        for header, seq in entries:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
        for header, seq in entries:
            acc = header.split()[0]
            h = int(hashlib.md5(acc.encode()).hexdigest()[:8], 16)
            dec = shuffle_protein(seq, seed ^ h, keep_cleavage_sites)
            out.write(f">{decoy_tag}{header}\n")
            for i in range(0, len(dec), 60):
                out.write(dec[i:i + 60] + "\n")
            n += 1
    return n, n


def split_fasta(fasta_in: str | Path, out_dir: str | Path, n_chunks: int) -> list[Path]:
    """Round-robin split of a FASTA into n_chunks files (for memory-bounded chunked searches)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    handles = [open(out_dir / f"chunk_{i:03d}.fasta", "w") for i in range(n_chunks)]
    try:
        for k, (header, seq) in enumerate(read_fasta(fasta_in)):
            fh = handles[k % n_chunks]
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    finally:
        for fh in handles:
            fh.close()
    return [out_dir / f"chunk_{i:03d}.fasta" for i in range(n_chunks)]
