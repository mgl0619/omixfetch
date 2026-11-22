"""Self-contained HTML summary of a download.

One file, no external assets or scripts, so it opens from a cluster scratch
directory as easily as from a laptop. Written for a dry run (showing the
plan) and again when a run finishes (showing outcomes).
"""

from __future__ import annotations

import html
import time
from pathlib import Path

from . import __version__

TREE = """run
|-- no public files ........................ skip (controlled access)
|-- submitted as BAM/CRAM  (header read with a 4 MB range request)
|   |-- 10x_bam_to_fastq lines ............. AWS BAM -> bamtofastq
|   |-- single-cell, no 10x lines .......... no faithful route
|   |-- aligner kept unmapped reads ........ AWS BAM -> samtools fastq
|   `-- dropped unmapped / unknown ......... no faithful route
`-- submitted as FASTQ
    `-- ENA FASTQ or AWS .sra -> fasterq-dump, fastest by measured throughput"""

CSS = """
body{font:14px/1.45 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:2rem auto;
     max-width:1200px;color:#1d2327;padding:0 1rem}
h1{font-size:1.5rem;margin:0 0 .2rem}h2{font-size:1.1rem;margin:2rem 0 .6rem;
   border-bottom:1px solid #dcdcde;padding-bottom:.3rem}
.sub{color:#646970;margin:0 0 1.2rem}
.cards{display:flex;gap:.8rem;flex-wrap:wrap}
.card{border:1px solid #dcdcde;border-radius:6px;padding:.6rem 1rem;min-width:120px}
.card b{display:block;font-size:1.3rem}.card span{color:#646970;font-size:.85rem}
table{border-collapse:collapse;width:100%;font-size:.88rem}
th,td{border-bottom:1px solid #eee;padding:.35rem .5rem;text-align:left;vertical-align:top}
th{background:#f6f7f7;font-weight:600}td.num{text-align:right;font-variant-numeric:tabular-nums}
.done{color:#007017;font-weight:600}.failed{color:#b32d2e;font-weight:600}
.skipped,.pending{color:#996800;font-weight:600}
.small{color:#646970;font-size:.8rem}
pre{background:#f6f7f7;padding:.8rem 1rem;border-radius:6px;font-size:.82rem;overflow-x:auto}
code{font-size:.85em}
"""


def _e(value) -> str:
    return html.escape("" if value is None else str(value))


def _bytes(n) -> str:
    if not n:
        return ""
    size = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024
    return ""


def _duration(seconds) -> str:
    if seconds is None:
        return ""
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h {m:02d}m" if h else (f"{m}m {s:02d}s" if m else f"{s}s")


def _spots(expected, written) -> str:
    if not expected and written is None:
        return ""
    if expected and written is not None:
        mark = "&#10003;" if expected == written else "&#10007;"
        return f"{expected:,} / {written:,} {mark}"
    if expected:
        return f"{expected:,} / &ndash;"
    return f"&ndash; / {written:,}"


def write_report(
    path: Path,
    accession: str,
    plans: list,
    store,
    *,
    title: str | None = None,
    rates: dict | None = None,
    tools: dict | None = None,
    dry_run: bool = False,
    wall_seconds: float | None = None,
) -> Path:
    """Write the HTML report and return its path."""
    states = {p.run: store.get(p.run) for p in plans}
    counts = {k: 0 for k in ("done", "failed", "skipped", "pending")}
    for s in states.values():
        counts[s.status if s.status in counts else "pending"] += 1
    moved = sum(s.bytes or 0 for s in states.values() if s.status == "done")
    busy = sum((s.finished - s.started) for s in states.values()
               if s.status == "done" and s.started and s.finished)

    rows = []
    for plan in plans:
        s = states[plan.run]
        d = plan.decision
        chosen = s.route if s.status == "done" else (d.order[0] if d and d.order else "")
        est = ""
        if d and chosen in d.estimates:
            est = f"{d.estimates[chosen]:.1f} h"
        excluded = "<br>".join(
            f"<code>{_e(r)}</code>: {_e(why)}" for r, why in
            (d.excluded.items() if d else (s.excluded or {}).items())
            if plan.routes.get(r)
        )
        took = (s.finished - s.started) if (s.started and s.finished) else None
        rate = f"{(s.bytes or 0) / took / 1024**2:.1f}" if took else ""
        files = "<br>".join(_e(Path(f).name) for f in s.fastq)
        rows.append(
            f"<tr><td>{_e(plan.run)}</td><td>{_e(plan.sample)}</td>"
            f"<td>{_e(plan.submitted)}</td><td>{_e(d.branch if d else s.branch)}</td>"
            f"<td><code>{_e(chosen)}</code></td><td class='small'>{excluded}</td>"
            f"<td class='num'>{_e(est)}</td>"
            f"<td class='{_e(s.status)}'>{_e(s.status)}</td>"
            f"<td class='num'>{_spots(plan.expected_spots, s.spots_written)}</td>"
            f"<td class='num'>{_bytes(s.bytes)}</td><td class='num'>{_duration(took)}</td>"
            f"<td class='num'>{rate}</td><td class='small'>{files}</td></tr>"
        )

    failures = [
        f"<tr><td>{_e(run)}</td><td class='{_e(s.status)}'>{_e(s.status)}</td>"
        f"<td>{_e(', '.join(s.attempts))}</td><td>{_e(s.error)}</td></tr>"
        for run, s in states.items() if s.status in ("failed", "skipped")
    ]
    rate_rows = "".join(
        f"<tr><td>{_e(k.upper())}</td><td class='num'>{v:.1f} MiB/s</td></tr>"
        for k, v in sorted((rates or {}).items())
    ) or "<tr><td colspan=2 class='small'>not measured</td></tr>"
    tool_rows = "".join(
        f"<tr><td>{_e(k)}</td><td>{_e(v or 'not found')}</td></tr>"
        for k, v in (tools or {}).items()
    )

    head = "<tr><th>Run</th><th>Status</th><th>Routes tried</th><th>Error</th></tr>"
    problems = (
        f"<table><thead>{head}</thead><tbody>{''.join(failures)}</tbody></table>"
        if failures else "<p>None.</p>"
    )
    kind = "Plan (dry run)" if dry_run else "Download report"
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>{_e(accession)} &ndash; omixfetch {kind.lower()}</title><style>{CSS}</style></head><body>
<h1>{_e(accession)} &ndash; {kind}</h1>
<p class="sub">{_e(title)}<br>omixfetch {_e(__version__)} &middot;
generated {_e(time.strftime('%Y-%m-%d %H:%M:%S %Z'))}</p>

<div class="cards">
<div class="card"><b>{len(plans)}</b><span>runs</span></div>
<div class="card"><b class="done">{counts['done']}</b><span>complete</span></div>
<div class="card"><b class="failed">{counts['failed']}</b><span>failed</span></div>
<div class="card"><b class="skipped">{counts['skipped']}</b><span>skipped</span></div>
<div class="card"><b class="pending">{counts['pending']}</b><span>pending</span></div>
<div class="card"><b>{_bytes(moved) or '0'}</b><span>transferred</span></div>
<div class="card"><b>{_duration(wall_seconds) or '&ndash;'}</b><span>wall time</span></div>
<div class="card"><b>{_duration(busy) or '&ndash;'}</b><span>run time (summed)</span></div>
</div>

<h2>Per-run decisions and outcomes</h2>
<table><thead><tr><th>Run</th><th>Sample</th><th>Submitted</th><th>Decision</th>
<th>Route</th><th>Not used, and why</th><th>Est.</th><th>Status</th>
<th>Spots expected / written</th><th>Size</th><th>Time</th><th>MiB/s</th><th>FASTQ</th>
</tr></thead><tbody>{''.join(rows)}</tbody></table>

<h2>Problems</h2>
{problems}

<h2>Decision tree</h2>
<pre>{_e(TREE)}</pre>
<p class="small">Routes derived from a submitted BAM (the .sra and ENA's FASTQ) carry
whatever the BAM lost; for 10x data that is the barcode/UMI read. Spot counts are
checked against NCBI RunInfo where the converter reports them.</p>

<h2>Measured throughput</h2>
<table><tbody>{rate_rows}</tbody></table>

<h2>Tools</h2>
<table><tbody>{tool_rows}</tbody></table>
</body></html>
"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(doc)
    return path
