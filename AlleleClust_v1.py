#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AlleleClust (single-file: core + GUI)

Purpose:
- Cluster alleles by exact CDS nucleotide sequence identity (FASTA from BAKO).
- Use BAKO SUMMARY TSV as the ONLY source of genus/species/replicon annotations.
- Output exactly 3 TSV files with your requested headers and rep formatting.

Run GUI:
    python AlleleClust.py

Run CLI:
    python AlleleClust.py --fasta "msrA offline_COMPLETE_CDS_nt.fasta" \
                          --bako-tsv "SUMMARY_COMPLETE.tsv" \
                          --outdir "ALLELECLUST_OUT" \
                          --prefix "msrA_offline" \
                          --zero-pad 3
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import traceback
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ============================================================
# CORE
# ============================================================

@dataclass(frozen=True)
class Member:
    accession: str
    header: str
    seq: str
    nt_length: int
    genus: str
    species: str
    replicon: str


def ensure_dir(path: str) -> None:
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Return list of (header_without_>, sequence)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA not found: {path}")

    out: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(chunks).replace(" ", "").upper()
                    out.append((header, seq))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())

    if header is not None:
        seq = "".join(chunks).replace(" ", "").upper()
        out.append((header, seq))

    return out


def read_bako_tsv(path: str) -> Dict[str, Dict[str, str]]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"BAKO TSV not found: {path}")

    with open(path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("BAKO TSV missing header row.")

        required = {"accession", "genus", "species", "replicon_type"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"BAKO TSV missing required columns: {sorted(missing)}")

        meta: Dict[str, Dict[str, str]] = {}
        dup = 0
        for row in reader:
            acc = (row.get("accession") or "").strip()
            if not acc:
                continue
            if acc in meta:
                dup += 1
            meta[acc] = row

        if dup:
            print(f"WARNING: {dup} duplicate accessions in TSV; last row wins.")

        return meta


def _is_blank_or_unknown(s: str) -> bool:
    x = (s or "").strip()
    if x == "":
        return True
    lx = x.lower()
    return lx in {"unknown", "na", "nan", "none", "null", "-"}


def _norm_replicon(s: str) -> str:
    x = (s or "").strip().lower()
    if x in {"chromosome", "plasmid", "phage", "prophage", "unknown"}:
        return x
    if x == "":
        return "unknown"
    return x


def cluster_exact(members: List[Member], allele_prefix: str = "allele", zero_pad: int = 0) -> Dict[str, List[Member]]:
    """Cluster members by exact CDS string identity."""
    seq_to_allele: Dict[str, str] = {}
    allele_to_members: Dict[str, List[Member]] = {}
    allele_i = 0

    for m in members:
        if not m.seq:
            continue

        aid = seq_to_allele.get(m.seq)
        if aid is None:
            allele_i += 1
            num = f"{allele_i:0{zero_pad}d}" if zero_pad > 0 else str(allele_i)
            aid = f"{allele_prefix}_{num}"
            seq_to_allele[m.seq] = aid
            allele_to_members[aid] = []

        allele_to_members[aid].append(m)

    return allele_to_members


def _rep_label(genus: str, species: str, accession: str) -> str:
    g = (genus or "").strip()
    s = (species or "").strip()

    if not _is_blank_or_unknown(g) and not _is_blank_or_unknown(s):
        return f"{g} {s}-{accession}"
    if not _is_blank_or_unknown(g):
        return f"{g}-{accession}"
    return accession


def _choose_min_accession(group: List[Member]) -> str:
    return min(m.accession for m in group)


def _build_species_reps(members: List[Member]) -> List[Tuple[str, str]]:
    """One rep per unique known (genus,species). Returns (label, accession)."""
    buckets: Dict[Tuple[str, str], List[Member]] = {}
    for m in members:
        g = (m.genus or "").strip()
        s = (m.species or "").strip()
        if _is_blank_or_unknown(g) or _is_blank_or_unknown(s):
            continue
        buckets.setdefault((g, s), []).append(m)

    reps: List[Tuple[str, str]] = []
    for (g, s), group in sorted(buckets.items(), key=lambda x: (x[0][0].lower(), x[0][1].lower())):
        acc = _choose_min_accession(group)
        reps.append((_rep_label(g, s, acc), acc))
    return reps


def _build_genus_reps_for_unknown_species(members: List[Member]) -> List[Tuple[str, str]]:
    """One rep per genus that has blank/unknown species. Returns (label, accession)."""
    buckets: Dict[str, List[Member]] = {}
    for m in members:
        g = (m.genus or "").strip()
        s = (m.species or "").strip()
        if _is_blank_or_unknown(g):
            continue
        if _is_blank_or_unknown(s):
            buckets.setdefault(g, []).append(m)

    reps: List[Tuple[str, str]] = []
    for g, group in sorted(buckets.items(), key=lambda x: x[0].lower()):
        acc = _choose_min_accession(group)
        reps.append((_rep_label(g, "", acc), acc))
    return reps


def _join(items: Sequence[str]) -> str:
    return ";".join(items)


def write_tsv(path: str, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

# ================== FASTA OUTPUT HELPERS ====================
def _safe_name(s: str) -> str:
    bad = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']
    out = s
    for ch in bad:
        out = out.replace(ch, "_")
    return out.strip().replace(" ", "_")



def write_fasta(path: str, records: Sequence[Tuple[str, str]]) -> None:
    """Write FASTA given (header_without_>, seq) records."""
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


# ================== REQUIRED FASTA OUTPUTS ====================
def _append_allele_to_header(header: str, allele_id: str) -> str:
    """
    Replace the first token (accession) in a FASTA header with:
      accession__allele_id
    Keeps the rest of the header unchanged.
    """
    h = (header or "").strip()
    if not h:
        return h
    parts = h.split(maxsplit=1)
    acc = parts[0]
    rest = parts[1] if len(parts) == 2 else ""
    new_acc = f"{acc}__{allele_id}"
    return f"{new_acc} {rest}".strip()


def write_required_fasta_outputs(
    alleles: Dict[str, List[Member]],
    outdir: str,
    prefix: str,
    reps_by_allele: Dict[str, Sequence[str]],
    logger=None,
) -> None:
    """
    ALWAYS writes:
      1) <prefix>_ALL_WITH_ALLELE_ID.fasta
         All members, but accession token becomes accession__allele_id

      2) <prefix>_REPRESENTATIVES_ONE_PER_ALLELE.fasta
         Exactly one representative per allele (deterministic).
         If reps_by_allele[allele_id] has multiple reps, choose min(accession).
         If empty, choose min(member accession) for that allele.
    """
    def log(msg: str):
        if logger:
            logger(msg)
        else:
            print(msg)

    # 1) All entries with allele_id appended
    all_records: List[Tuple[str, str]] = []
    for allele_id, members in alleles.items():
        for m in members:
            new_header = _append_allele_to_header(m.header, allele_id)
            all_records.append((new_header, m.seq))

    all_fp = os.path.join(outdir, f"{prefix}_ALL_WITH_ALLELE_ID.fasta")
    write_fasta(all_fp, all_records)

    # 2) One representative per allele (deterministic)
    rep_records: List[Tuple[str, str]] = []
    for allele_id, members in alleles.items():
        if not members:
            continue

        wanted = list(reps_by_allele.get(allele_id, []))
        by_acc = {m.accession: m for m in members}

        chosen_member: Optional[Member] = None

        if wanted:
            chosen_acc = min(wanted)
            chosen_member = by_acc.get(chosen_acc)

        # Fallback: no reps or rep not present for some reason
        if chosen_member is None:
            chosen_acc = min(m.accession for m in members)
            chosen_member = by_acc[chosen_acc]

        rep_records.append((chosen_member.header, chosen_member.seq))

    rep_fp = os.path.join(outdir, f"{prefix}_REPRESENTATIVES_ONE_PER_ALLELE.fasta")
    write_fasta(rep_fp, rep_records)

    log("✅ Required FASTA outputs written:")
    log(f"  All members (allele-tagged accessions): {all_fp}")
    log(f"  One rep per allele: {rep_fp}")


def write_optional_fasta_outputs(
    alleles: Dict[str, List[Member]],
    outdir: str,
    prefix: str,
    reps_by_allele: Dict[str, Sequence[str]],
    logger=None,
) -> None:
    """Write optional FASTA folders:
    - <prefix>_ALLELE_FASTAS: one FASTA per allele containing all member CDS
    - <prefix>_REP_FASTAS: one FASTA per allele containing representative CDS (may be multiple)
    reps_by_allele maps allele_id -> iterable of accession strings to include.
    """
    def log(msg: str):
        if logger:
            logger(msg)
        else:
            print(msg)

    allele_dir = os.path.join(outdir, f"{prefix}_ALLELE_FASTAS")
    rep_dir = os.path.join(outdir, f"{prefix}_REP_FASTAS")
    ensure_dir(allele_dir)
    ensure_dir(rep_dir)

    for allele_id, members in alleles.items():
        if not members:
            continue

        allele_fp = os.path.join(allele_dir, f"{_safe_name(allele_id)}.fasta")
        group_records = [(m.header, m.seq) for m in members]
        write_fasta(allele_fp, group_records)

        wanted = set(reps_by_allele.get(allele_id, []))
        if wanted:
            by_acc = {m.accession: m for m in members}
            rep_records = []
            for acc in sorted(wanted):
                m = by_acc.get(acc)
                if m:
                    rep_records.append((m.header, m.seq))
            if rep_records:
                rep_fp = os.path.join(rep_dir, f"{_safe_name(allele_id)}_representatives.fasta")
                write_fasta(rep_fp, rep_records)

    log("✅ Optional FASTA outputs written:")
    log(f"  Allele groups: {allele_dir}")
    log(f"  Representatives: {rep_dir}")


def run_alleleclust(
    fasta: str,
    bako_tsv: str,
    outdir: str,
    prefix: str = "alleles",
    allele_prefix: str = "allele",
    zero_pad: int = 0,
    write_fastas: bool = False,
    logger=None,  # optional callable for GUI logging
) -> None:
    """Main entry point used by both CLI and GUI."""
    def log(msg: str):
        if logger:
            logger(msg)
        else:
            print(msg)

    ensure_dir(outdir)
    fasta_records = read_fasta(fasta)
    meta = read_bako_tsv(bako_tsv)

    members: List[Member] = []
    fasta_accs: List[str] = []

    for header, seq in fasta_records:
        if not header:
            continue
        acc = header.split()[0].strip()
        fasta_accs.append(acc)

        m = meta.get(acc, {})
        genus = (m.get("genus") or "").strip()
        species = (m.get("species") or "").strip()
        replicon = _norm_replicon(m.get("replicon_type") or "unknown")

        members.append(Member(acc, header, seq, len(seq), genus, species, replicon))

    if not members:
        raise RuntimeError("No FASTA records found.")

    # Coverage warning (prevents silent mismatches)
    fs, ts = set(fasta_accs), set(meta.keys())
    matched = fs & ts
    rate = (len(matched) / len(fs) * 100.0) if fs else 0.0
    if rate < 99.0:
        log(f"WARNING: match rate FASTA↔TSV {rate:.2f}% ({len(matched)}/{len(fs)})")
        log(f"  FASTA not in TSV (first 25): {sorted(list(fs-ts))[:25]}")
        log(f"  TSV not in FASTA (first 25): {sorted(list(ts-fs))[:25]}")

    alleles = cluster_exact(members, allele_prefix=allele_prefix, zero_pad=zero_pad)

    # Output 2: members
    member_rows: List[Dict[str, str]] = []
    for aid, ms in alleles.items():
        for m in ms:
            member_rows.append({
                "accession": m.accession,
                "allele_id": aid,
                "header": m.header,
                "nt_length": str(m.nt_length),
                "genus": m.genus,
                "species": m.species,
                "replicon": m.replicon,
            })

    # Output 1: allele description
    desc_rows: List[Dict[str, str]] = []
    for aid, ms in alleles.items():
        if not ms:
            continue
        nt_len = ms[0].nt_length
        aa_len = nt_len // 3
        accs = sorted({m.accession for m in ms})
        example_header = ms[0].header

        genera = sorted({(m.genus or "").strip() for m in ms if not _is_blank_or_unknown(m.genus)}, key=lambda x: x.lower())
        species_pairs = sorted({
            f"{(m.genus or '').strip()} {(m.species or '').strip()}".strip()
            for m in ms
            if (not _is_blank_or_unknown(m.genus)) and (not _is_blank_or_unknown(m.species))
        }, key=lambda x: x.lower())
        replicons = sorted({_norm_replicon(m.replicon) for m in ms if (m.replicon or "").strip() != ""})

        desc_rows.append({
            "allele_id": aid,
            "nt_length": str(nt_len),
            "aa_length": str(aa_len),
            "n_members": str(len(ms)),
            "accessions": _join(accs),
            "example_header": example_header,
            "present_genus": _join(genera),
            "present_species": _join(species_pairs),
            "present_replicon": _join(replicons),
            "genus_count": str(len(genera)),
            "species_count": str(len(species_pairs)),
        })

    # Output 3: representatives
    reps_rows: List[Dict[str, str]] = []
    reps_accessions_by_allele: Dict[str, List[str]] = {}
    for aid, ms in alleles.items():
        if not ms:
            continue

        genus_presence = sorted({(m.genus or "").strip() for m in ms if not _is_blank_or_unknown(m.genus)}, key=lambda x: x.lower())
        species_presence = sorted({
            f"{(m.genus or '').strip()} {(m.species or '').strip()}".strip()
            for m in ms
            if (not _is_blank_or_unknown(m.genus)) and (not _is_blank_or_unknown(m.species))
        }, key=lambda x: x.lower())

        allele_genus = genus_presence[0] if len(genus_presence) == 1 else ("mixed" if len(genus_presence) > 1 else "")
        repl_set = sorted({_norm_replicon(m.replicon) for m in ms if (m.replicon or "").strip() != ""})
        allele_type = "+".join(repl_set) if repl_set else "unknown"

        species_reps_t = _build_species_reps(ms)
        genus_reps_t = _build_genus_reps_for_unknown_species(ms)
        species_reps = [x[0] for x in species_reps_t]
        genus_reps = [x[0] for x in genus_reps_t]
        rep_accs_all = {x[1] for x in (species_reps_t + genus_reps_t)}

        def reps_for_repl(target: str) -> Tuple[str, set]:
            subset = [m for m in ms if _norm_replicon(m.replicon) == target]
            if not subset:
                return "", set()
            s_t = _build_species_reps(subset)
            g_t = _build_genus_reps_for_unknown_species(subset)
            labels = [x[0] for x in (s_t + g_t)]
            accs = {x[1] for x in (s_t + g_t)}
            return _join(labels), accs

        chromosome_rep, chr_accs = reps_for_repl("chromosome")
        plasmid_rep, pls_accs = reps_for_repl("plasmid")
        phage_rep, phg_accs = reps_for_repl("phage")
        rep_accs_all |= (chr_accs | pls_accs | phg_accs)

        reps_rows.append({
            "allele_id": aid,
            "allele_type": allele_type,
            "allele_genus": allele_genus,
            "genus_presence": _join(genus_presence),
            "species_presence": _join(species_presence),
            "genus_rep": _join(genus_reps),
            "species_rep": _join(species_reps),
            "species_count": str(len(species_presence)),
            "chromosome_rep": chromosome_rep,
            "plasmid_rep": plasmid_rep,
            "phage_rep": phage_rep,
            "member_count": str(len(ms)),
        })
        reps_accessions_by_allele[aid] = sorted(rep_accs_all)

    members_path = os.path.join(outdir, f"{prefix}_ALLELE_MEMBERS_FINAL.tsv")
    desc_path = os.path.join(outdir, f"{prefix}_ALLELES_DESCRIPTION_FINAL.tsv")
    reps_path = os.path.join(outdir, f"{prefix}_REPRESENTATIVES_FINAL.tsv")

    write_tsv(members_path,
              ["accession", "allele_id", "header", "nt_length", "genus", "species", "replicon"],
              member_rows)
    write_tsv(desc_path,
              ["allele_id", "nt_length", "aa_length", "n_members", "accessions", "example_header",
               "present_genus", "present_species", "present_replicon", "genus_count", "species_count"],
              desc_rows)
    write_tsv(reps_path,
              ["allele_id", "allele_type", "allele_genus", "genus_presence", "species_presence",
               "genus_rep", "species_rep", "species_count", "chromosome_rep", "plasmid_rep", "phage_rep", "member_count"],
              reps_rows)

    if write_fastas:
        write_optional_fasta_outputs(
            alleles=alleles,
            outdir=outdir,
            prefix=prefix,
            reps_by_allele=reps_accessions_by_allele,
            logger=logger,
        )

    # ALWAYS write these two FASTA outputs
    write_required_fasta_outputs(
        alleles=alleles,
        outdir=outdir,
        prefix=prefix,
        reps_by_allele=reps_accessions_by_allele,
        logger=logger,
    )

    log("✅ Run complete.")
    log(f"Output folder: {outdir}")


# ============================================================
# GUI (PyQt6) — loaded only when needed
# ============================================================

def launch_gui() -> None:
    try:
        from PyQt6.QtCore import QThread, pyqtSignal
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
            QLabel, QLineEdit, QPushButton, QFileDialog, QMessageBox,
            QSpinBox, QPlainTextEdit, QGroupBox, QCheckBox
        )
    except Exception as e:
        print("ERROR: PyQt6 is not installed or not available in this environment.")
        print("Install with: pip install PyQt6")
        raise

    class Worker(QThread):
        log = pyqtSignal(str)
        done = pyqtSignal()
        failed = pyqtSignal(str)

        def __init__(self, fasta, bako_tsv, outdir, prefix, allele_prefix, zero_pad, write_fastas):
            super().__init__()
            self.fasta = fasta
            self.bako_tsv = bako_tsv
            self.outdir = outdir
            self.prefix = prefix
            self.allele_prefix = allele_prefix
            self.zero_pad = zero_pad
            self.write_fastas = write_fastas

        def run(self):
            try:
                self.log.emit("Starting AlleleClust…")
                self.log.emit(f"FASTA: {self.fasta}")
                self.log.emit(f"BAKO TSV: {self.bako_tsv}")
                self.log.emit(f"OUTDIR: {self.outdir}")
                self.log.emit(f"prefix={self.prefix} allele_prefix={self.allele_prefix} zero_pad={self.zero_pad}")
                self.log.emit(f"write_fastas={self.write_fastas}")
                self.log.emit("-" * 70)

                run_alleleclust(
                    fasta=self.fasta,
                    bako_tsv=self.bako_tsv,
                    outdir=self.outdir,
                    prefix=self.prefix,
                    allele_prefix=self.allele_prefix,
                    zero_pad=self.zero_pad,
                    write_fastas=self.write_fastas,
                    logger=self.log.emit,
                )

                self.done.emit()
            except Exception:
                self.failed.emit(traceback.format_exc())

    class GUI(QWidget):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("AlleleClust (BAKO FASTA + BAKO TSV)")
            self.setMinimumWidth(860)

            layout = QVBoxLayout(self)

            # Inputs
            box = QGroupBox("Inputs")
            form = QFormLayout(box)

            self.fasta_path = QLineEdit()
            btn_fasta = QPushButton("Browse")
            btn_fasta.clicked.connect(self._browse_fasta)
            row = QHBoxLayout()
            row.addWidget(self.fasta_path)
            row.addWidget(btn_fasta)
            form.addRow(QLabel("FASTA (COMPLETE_CDS_nt.fasta):"), self._wrap(row))

            self.tsv_path = QLineEdit()
            btn_tsv = QPushButton("Browse")
            btn_tsv.clicked.connect(self._browse_tsv)
            row = QHBoxLayout()
            row.addWidget(self.tsv_path)
            row.addWidget(btn_tsv)
            form.addRow(QLabel("BAKO SUMMARY TSV:"), self._wrap(row))

            self.out_dir = QLineEdit()
            btn_out = QPushButton("Browse")
            btn_out.clicked.connect(self._browse_outdir)
            row = QHBoxLayout()
            row.addWidget(self.out_dir)
            row.addWidget(btn_out)
            form.addRow(QLabel("Output folder:"), self._wrap(row))

            layout.addWidget(box)

            # Options
            opt = QGroupBox("Options")
            opt_form = QFormLayout(opt)

            self.prefix = QLineEdit("msrA_offline")
            opt_form.addRow(QLabel("Output prefix:"), self.prefix)

            self.allele_prefix = QLineEdit("allele")
            opt_form.addRow(QLabel("Allele ID prefix:"), self.allele_prefix)

            self.zero_pad = QSpinBox()
            self.zero_pad.setRange(0, 9)
            self.zero_pad.setValue(0)
            opt_form.addRow(QLabel("Zero-pad allele numbers:"), self.zero_pad)

            # FASTA output checkbox
            self.write_fastas = QCheckBox("Write FASTA folders (allele groups + representatives)")
            self.write_fastas.setChecked(False)
            opt_form.addRow(self.write_fastas)

            layout.addWidget(opt)

            # Run
            self.btn_run = QPushButton("RUN")
            self.btn_run.clicked.connect(self._run)
            layout.addWidget(self.btn_run)

            # Log
            layout.addWidget(QLabel("Log:"))
            self.log = QPlainTextEdit()
            self.log.setReadOnly(True)
            self.log.setMinimumHeight(260)
            layout.addWidget(self.log)

            self.worker = None

        def _wrap(self, hbox):
            from PyQt6.QtWidgets import QWidget
            w = QWidget()
            w.setLayout(hbox)
            return w

        def _browse_fasta(self):
            from PyQt6.QtWidgets import QFileDialog
            path, _ = QFileDialog.getOpenFileName(self, "Select FASTA", "", "FASTA (*.fasta *.fa *.fna *.txt);;All files (*)")
            if path:
                self.fasta_path.setText(path)

        def _browse_tsv(self):
            from PyQt6.QtWidgets import QFileDialog
            path, _ = QFileDialog.getOpenFileName(self, "Select BAKO SUMMARY TSV", "", "TSV (*.tsv);;All files (*)")
            if path:
                self.tsv_path.setText(path)

        def _browse_outdir(self):
            from PyQt6.QtWidgets import QFileDialog
            path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
            if path:
                self.out_dir.setText(path)

        def _log(self, s: str):
            self.log.appendPlainText(s)

        def _validate(self) -> bool:
            from PyQt6.QtWidgets import QMessageBox
            fasta = self.fasta_path.text().strip()
            tsv = self.tsv_path.text().strip()
            outdir = self.out_dir.text().strip()

            if not fasta or not os.path.exists(fasta):
                QMessageBox.warning(self, "Missing file", "Please choose a valid FASTA file.")
                return False
            if not tsv or not os.path.exists(tsv):
                QMessageBox.warning(self, "Missing file", "Please choose a valid BAKO SUMMARY TSV.")
                return False
            if not outdir:
                QMessageBox.warning(self, "Missing folder", "Please choose an output folder.")
                return False
            try:
                os.makedirs(outdir, exist_ok=True)
            except Exception as e:
                QMessageBox.critical(self, "Output folder error", str(e))
                return False

            if not self.prefix.text().strip():
                QMessageBox.warning(self, "Missing prefix", "Please set an output prefix.")
                return False
            if not self.allele_prefix.text().strip():
                QMessageBox.warning(self, "Missing allele prefix", "Please set an allele ID prefix.")
                return False
            return True

        def _run(self):
            from PyQt6.QtWidgets import QMessageBox
            if not self._validate():
                return

            fasta = self.fasta_path.text().strip()
            tsv = self.tsv_path.text().strip()
            outdir = self.out_dir.text().strip()
            prefix = self.prefix.text().strip()
            allele_prefix = self.allele_prefix.text().strip()
            zero_pad = int(self.zero_pad.value())
            write_fastas = bool(self.write_fastas.isChecked())

            self.btn_run.setEnabled(False)
            self._log("")
            self._log("=" * 70)

            self.worker = Worker(fasta, tsv, outdir, prefix, allele_prefix, zero_pad, write_fastas)
            self.worker.log.connect(self._log)
            self.worker.done.connect(lambda: self._on_done())
            self.worker.failed.connect(lambda tb: self._on_failed(tb))
            self.worker.start()

        def _on_done(self):
            from PyQt6.QtWidgets import QMessageBox
            self.btn_run.setEnabled(True)
            QMessageBox.information(self, "Done", "Run complete. Check the output folder.")

        def _on_failed(self, tb: str):
            from PyQt6.QtWidgets import QMessageBox
            self.btn_run.setEnabled(True)
            self._log("❌ RUN FAILED")
            self._log(tb)
            QMessageBox.critical(self, "Run failed", "Run failed. See log for traceback.")

    app = QApplication(sys.argv)
    w = GUI()
    w.show()
    sys.exit(app.exec())


# ============================================================
# ENTRYPOINT
# ============================================================

def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(add_help=True)

    # If no CLI args, we launch GUI.
    p.add_argument("--fasta", help="BAKO COMPLETE CDS nt FASTA")
    p.add_argument("--bako-tsv", help="BAKO SUMMARY TSV")
    p.add_argument("--outdir", help="Output folder")
    p.add_argument("--prefix", default="alleles", help="Output filename prefix")
    p.add_argument("--allele-prefix", default="allele", help="Allele ID prefix")
    p.add_argument("--zero-pad", type=int, default=0, help="Zero pad allele numbers, e.g. 3 -> allele_001")
    p.add_argument("--write-fastas", action="store_true", help="Write per-allele and representative FASTA folders")
    p.add_argument("--gui", action="store_true", help="Launch GUI explicitly")
    return p.parse_args(argv)


def main() -> None:
    args = parse_args(sys.argv[1:])

    # GUI if requested OR if user provided no required CLI inputs
    if args.gui or (not args.fasta and not args.bako_tsv and not args.outdir):
        launch_gui()
        return

    # CLI requires these three
    if not args.fasta or not args.bako_tsv or not args.outdir:
        print("ERROR: For CLI mode you must provide --fasta, --bako-tsv, and --outdir.")
        print("Tip: run without arguments to open the GUI.")
        sys.exit(2)

    run_alleleclust(
        fasta=args.fasta,
        bako_tsv=args.bako_tsv,
        outdir=args.outdir,
        prefix=args.prefix,
        allele_prefix=args.allele_prefix,
        zero_pad=args.zero_pad,
        write_fastas=args.write_fastas,
        logger=None,
    )


if __name__ == "__main__":
    main()
