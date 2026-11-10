#!/usr/bin/env python3
"""Assemble a read-accounting table.

Inputs: --sheet step1_samplesheet.tsv, --seqkit directory of <sample_id>.seqkit.tsv,
--logdir host-filter logs (optional; adds per-host removal columns), --out TSV.

Counts are total reads (R1+R2): raw and post-fastp from fastp JSON, post-host from
seqkit. Per-host removal is read from the "processed N reads" lines that each
samtools fastq call writes to stderr (cow, GRCh38, T2T-Y, in that order). For
samples with several job ids, the log whose last value equals post-host is used.
"""
import argparse, json, os, sys, csv, re, glob

_PROC_RE = re.compile(r"bam2fq_mainloop\] processed (\d+) reads")

def build_hostlog_index(logdir):
    """sample_id -> [jobid, ...] from logs/*.out 'Sample name:' lines."""
    idx = {}
    for outp in glob.glob(os.path.join(logdir, "*.out")):
        job = os.path.basename(outp)[:-4]
        try:
            with open(outp, errors="replace") as fh:
                for line in fh:
                    if line.startswith("Sample name:"):
                        idx.setdefault(line.split("Sample name:", 1)[1].strip(), []).append(job)
                        break
        except Exception:
            pass
    return idx

def _parse_err_processed(errp):
    vals = []
    try:
        with open(errp, errors="replace") as fh:
            for line in fh:
                m = _PROC_RE.search(line)
                if m:
                    vals.append(int(m.group(1)))
    except Exception:
        pass
    return vals

def host_breakdown(jobids, logdir, postfastp, posthost):
    """Return (cow_rm, grch38_rm, t2ty_rm, human_rm, mismatch_bool) or None.
    Picks the .err whose last 'processed' == posthost (handles re-runs)."""
    if postfastp is None or posthost is None:
        return None
    best = None
    for job in jobids:
        vals = _parse_err_processed(os.path.join(logdir, job + ".err"))
        if len(vals) >= 3 and vals[-1] == posthost:
            best = vals[-3:]
            break
    if best is None:
        return None
    cow_out, grch_out, t2t_out = best
    cow_rm = postfastp - cow_out
    grch_rm = cow_out - grch_out
    t2ty_rm = grch_out - t2t_out
    human_rm = grch_rm + t2ty_rm
    # telescopes exactly unless a pair-filtered re-run dropped a few reads
    mismatch = (cow_rm + human_rm) != (postfastp - posthost)
    return cow_rm, grch_rm, t2ty_rm, human_rm, mismatch

def read_fastp(path):
    try:
        with open(path) as fh:
            d = json.load(fh)
        s = d["summary"]
        return (int(s["before_filtering"]["total_reads"]),
                int(s["after_filtering"]["total_reads"]),
                int(s["before_filtering"].get("total_bases", 0)),
                int(s["after_filtering"].get("total_bases", 0)))
    except Exception as e:
        sys.stderr.write(f"WARN: could not parse fastp json {path}: {e}\n")
        return (None, None, None, None)

def read_seqkit(path):
    """Sum num_seqs and sum_len across rows; capture avg/max len from R1 row."""
    try:
        reads = bases = 0
        avg_len = max_len = None
        with open(path) as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for i, row in enumerate(r):
                reads += int(float(row["num_seqs"]))
                bases += int(float(row["sum_len"]))
                if i == 0:
                    avg_len = row.get("avg_len")
                    max_len = row.get("max_len")
        return reads, bases, avg_len, max_len
    except Exception as e:
        sys.stderr.write(f"WARN: could not parse seqkit tsv {path}: {e}\n")
        return (None, None, None, None)

def pct(num, den):
    return round(100.0 * num / den, 4) if (den and den > 0) else ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--seqkit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--logdir", default=None,
                    help="host_filter/logs dir (.out/.err); enables per-host removal columns")
    a = ap.parse_args()

    hostlog_idx = build_hostlog_index(a.logdir) if a.logdir else {}

    cols = ["sample_id","pnum","field_name","compartment","gdna_conc","pcr_conc",
            "raw_reads","postfastp_reads","posthost_reads",
            "pct_removed_fastp","host_reads_removed","pct_host_removed","pct_surviving",
            "cow_reads_removed","pct_cow_removed",
            "human_grch38_reads_removed","human_t2ty_reads_removed",
            "human_reads_removed","pct_human_removed",
            "posthost_bases","posthost_avg_len","posthost_max_len","flags"]
    rows = []
    with open(a.sheet) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("sample_id\t"):
                continue
            f = line.split("\t")
            # tolerate short control rows
            sid, pnum, field, comp = (f+[""]*4)[:4]
            fastp_json = f[6] if len(f) > 6 else ""
            gdna = f[7] if len(f) > 7 else ""
            pcr  = f[8] if len(f) > 8 else ""

            raw = pf = None
            if fastp_json and os.path.exists(fastp_json):
                raw, pf, _, _ = read_fastp(fastp_json)

            sk = os.path.join(a.seqkit, f"{sid}.seqkit.tsv")
            ph = phbases = avg = mx = None
            if os.path.exists(sk):
                ph, phbases, avg, mx = read_seqkit(sk)

            flags = []
            if raw is None: flags.append("no_fastp_json")
            if ph  is None: flags.append("no_seqkit")
            if (pf is not None and ph is not None and ph > pf):
                flags.append("posthost>postfastp?")  # sense check:shouldn't happen

            # per-host breakdown (only meaningful for samples with surviving reads)
            cow_rm = grch_rm = t2ty_rm = human_rm = ""
            pct_cow = pct_human = ""
            if a.logdir and pf is not None and ph is not None and ph > 0:
                hb = host_breakdown(hostlog_idx.get(sid, []), a.logdir, pf, ph)
                if hb is None:
                    flags.append("no_hostlog")
                else:
                    cow_rm, grch_rm, t2ty_rm, human_rm, mismatch = hb
                    pct_cow = pct(cow_rm, pf)
                    pct_human = pct(human_rm, pf)
                    if mismatch:
                        flags.append("host_split_approx")  # pair-filtered re-run

            rows.append([
                sid, pnum, field, comp, gdna, pcr,
                raw if raw is not None else "",
                pf  if pf  is not None else "",
                ph  if ph  is not None else "",
                pct((raw-pf), raw) if (raw is not None and pf is not None) else "",
                (pf-ph) if (pf is not None and ph is not None) else "",
                pct((pf-ph), pf) if (pf is not None and ph is not None) else "",
                pct(ph, raw) if (raw is not None and ph is not None) else "",
                cow_rm, pct_cow,
                grch_rm, t2ty_rm,
                human_rm, pct_human,
                phbases if phbases is not None else "",
                avg if avg is not None else "",
                mx if mx is not None else "",
                ";".join(flags),
            ])

    with open(a.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)

    # brief stderr summary by compartment
    from collections import defaultdict
    agg = defaultdict(lambda: [0,0,0,0])  # n, sum_raw, sum_posthost, n_flagged
    for r in rows:
        comp = r[3] or "unknown"
        agg[comp][0] += 1
        if r[6] != "": agg[comp][1] += int(r[6])
        if r[8] != "": agg[comp][2] += int(r[8])
        if r[-1] != "": agg[comp][3] += 1   # flags is always the last column
    sys.stderr.write("\n== Accounting summary by compartment ==\n")
    sys.stderr.write(f"{'compartment':<16}{'n':>5}{'sum_raw':>16}{'sum_posthost':>16}{'flagged':>9}\n")
    for comp in sorted(agg):
        n, sr, sp, nf = agg[comp]
        sys.stderr.write(f"{comp:<16}{n:>5}{sr:>16,}{sp:>16,}{nf:>9}\n")
    sys.stderr.write(f"\nWrote {len(rows)} rows -> {a.out}\n")

if __name__ == "__main__":
    main()
