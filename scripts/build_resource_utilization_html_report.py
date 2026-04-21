#!/usr/bin/env python3
"""Build an HTML resource utilization report – stage-first layout with sortable tables."""
import json, os
from collections import defaultdict
from datetime import datetime
import re as _re
import html as _html

ROOT = "/Users/cmarkello/work/pacbio_workflows/hifisomatic_out/hifisomatic_tn_chr22_nogpu/20260415_180705_hifisomatic"
OUT_HTML = os.path.join(ROOT, "resource_utilization_report.html")

recs = json.load(open("/tmp/resources_run2.json"))

# ---------------------------------------------------------------------------
# Stage / task classification
# ---------------------------------------------------------------------------
STAGE_MAP = {
    "align_tumor":            ("Alignment",            "align_single_bam (pbmm2 → tumor)"),
    "align_normal":           ("Alignment",            "align_single_bam (pbmm2 → normal)"),
    "MosdepthTumor":          ("Alignment",            "mosdepth (tumor)"),
    "MosdepthNormal":         ("Alignment",            "mosdepth (normal)"),
    "bamstatsTumor":          ("Alignment",            "seqkit_bamstats (tumor)"),
    "bamstatsNormal":         ("Alignment",            "seqkit_bamstats (normal)"),
    "summarize_tumor_RL":     ("Alignment",            "summarize_seqkit_alignment (tumor)"),
    "summarize_normal_RL":    ("Alignment",            "summarize_seqkit_alignment (normal)"),
    "split_contigs":          ("Alignment",            "split_contigs"),
    "run_deepsomatic":        ("Small-variant calling","deepsomatic_paired (shards)"),
    "gather_deepsomatic":     ("Small-variant calling","gather_deepsomatic"),
    "clair3_tumor":           ("Small-variant calling","clair3 (tumor germline)"),
    "clair3_normal":          ("Small-variant calling","clair3 (normal germline)"),
    "correct_clair3_tumor":   ("Small-variant calling","correct_vcf (tumor)"),
    "correct_clair3_normal":  ("Small-variant calling","correct_vcf (normal)"),
    "phaseTumorBam":          ("Phasing",              "hiphase_with_somatic (tumor)"),
    "phaseNormalBam":         ("Phasing",              "hiphase (normal)"),
    "normSomatic":            ("Annotation",           "bcftools_norm (somatic)"),
    "annotateSomatic":        ("Annotation",           "vep_annotate (somatic)"),
    "tmb_estimate":           ("Annotation",           "tmb_estimate"),
    "prioritizeSomatic":      ("Annotation",           "prioritize_small_variants"),
    "mutationalpattern":      ("Annotation",           "mutationalpattern"),
    "pileup_tumor":           ("Methylation",          "cpg_pileup (tumor)"),
    "pileup_normal":          ("Methylation",          "cpg_pileup (normal)"),
    "DSS_DMR":                ("Methylation",          "DSS_DMR"),
    "annotate_DMR":           ("Methylation",          "annotate_DMR"),
    "prioritize_dmr_intogen": ("Methylation",          "prioritize_dmr_intogen"),
    "msi_score":              ("Biomarkers",           "owl_msi"),
    "phased_severus":         ("Structural variants",  "Severus_sv (phased)"),
    "unphased_severus":       ("Structural variants",  "Severus_sv (unphased)"),
    "wakhan":                 ("Structural variants",  "wakhan"),
    "SAVANA_sv":              ("Structural variants",  "SAVANA_sv"),
    "tabixSeverus":           ("Structural variants",  "tabix_vcf (Severus)"),
    "filter_Severus":         ("Structural variants",  "svpack_filter_annotated"),
    "recover_mate_bnd":       ("Structural variants",  "recover_mate_bnd"),
    "circos_BND":             ("Structural variants",  "circos_BND"),
    "annotateSeverus":        ("Structural variants",  "annotsv"),
    "prioritize_Severus":     ("Structural variants",  "prioritize_sv_intogen"),
    "cnvkit_tumor":           ("CNV / clonality",      "cnvkit_tumor"),
    "cnvkit_recall":          ("CNV / clonality",      "cnvkit_recall"),
    "mergeGermline":          ("CNV / clonality",      "merge_germline"),
    "Amber":                  ("CNV / clonality",      "Amber"),
    "Cobalt":                 ("CNV / clonality",      "Cobalt"),
    "purple_nosomatic":       ("CNV / clonality",      "Purple (nosomatic)"),
    "purple_withsomatic":     ("CNV / clonality",      "Purple (withsomatic)"),
    "chord_hrd":              ("CNV / clonality",      "chord_hrd"),
    "report_sample":          ("Reporting",            "report_sample"),
}

STAGE_ORDER = [
    "Alignment", "Small-variant calling", "Phasing", "Annotation",
    "Methylation", "Biomarkers", "Structural variants", "CNV / clonality",
    "Reporting", "Other",
]

STAGE_COLOR = {
    "Alignment":            "#4c78a8",
    "Small-variant calling":"#f58518",
    "Phasing":              "#e45756",
    "Annotation":           "#72b7b2",
    "Methylation":          "#54a24b",
    "Biomarkers":           "#b279a2",
    "Structural variants":  "#ff9da6",
    "CNV / clonality":      "#9d755d",
    "Reporting":            "#bab0ac",
    "Other":                "#aec7e8",
}

def classify(rec):
    call = rec["call"]
    if "deepsomatic_paired" in call:
        return STAGE_MAP["run_deepsomatic"]
    if "owl_score" in call or "owl_profile" in call:
        return ("Biomarkers", "owl_msi (" + call.rsplit("call-", 1)[-1].split("/")[0] + ")")
    if "align_single_bam" in call:
        if "align_tumor" in call:
            return STAGE_MAP["align_tumor"]
        if "align_normal" in call:
            return STAGE_MAP["align_normal"]
    name = call.split("/")[-1]
    if name.startswith("call-"):
        name = name[5:]
    name_plain = name.rsplit("-", 1)[0] if name.rsplit("-", 1)[-1].isdigit() else name
    for key in (name_plain, name):
        if key in STAGE_MAP:
            return STAGE_MAP[key]
    return ("Other", name)

# ---------------------------------------------------------------------------
# Scaling assumptions (chr22 → WGS)
# ---------------------------------------------------------------------------
CHR22_BP   = 50_818_468
WGS_BP     = 3_088_269_832
GENOME_SCALE = WGS_BP / CHR22_BP   # ≈ 60.77

def scale_for_rec(rec):
    call, task = rec["call"], rec["task_name"]
    g = GENOME_SCALE
    if "align_single_bam" in call:            return (g, 1.1, g * 0.9)
    if task == "mosdepth":                    return (g, 1.1, g * 0.9)
    if task in ("seqkit_bamstats", "summarize_seqkit_alignment"):
                                              return (g, 1.1, g)
    if task == "split_contigs":               return (1, 1, 1)
    if "deepsomatic_paired" in call:          return (g, 1.2, 3.0)
    if "gather_deepsomatic" in call:          return (g, 1.2, g)
    if "clair3_tumor" in call or "clair3_normal" in call:
                                              return (g, 1.2, g * 0.9)
    if "correct_clair3" in call:              return (g, 1.0, g)
    if task in ("hiphase", "hiphase_with_somatic"):
                                              return (g, 1.2, g * 0.9)
    if task in ("bcftools_norm", "vep_annotate", "tmb_estimate",
                "prioritize_small_variants", "mutationalpattern"):
                                              return (g, 1.2, g)
    if task == "cpg_pileup":                  return (g, 1.2, g * 0.9)
    if task == "DSS_DMR":                     return (g, 1.5, g)
    if task in ("annotate_DMR", "prioritize_dmr_intogen"):
                                              return (g, 1.0, g * 0.8)
    if "owl_profile" in call or "owl_score" in call:
                                              return (g, 1.0, g)
    if task in ("Severus_sv", "SAVANA_sv", "wakhan"):
                                              return (g, 1.5, g * 0.9)
    if task in ("tabix_vcf", "svpack_filter_annotated",
                "recover_mate_bnd", "circos_BND"):
                                              return (g * 0.5, 1.0, g * 0.5)
    if task == "annotsv":                     return (g, 1.5, g)
    if task == "prioritize_sv_intogen":       return (g, 1.0, g)
    if task == "cnvkit_tumor":                return (g, 1.5, g * 0.9)
    if task in ("cnvkit_recall", "merge_germline"):
                                              return (g, 1.2, g)
    if task in ("Amber", "Cobalt"):           return (g, 1.3, g * 0.8)
    if task == "Purple":                      return (g, 1.3, g)
    if task == "chord_hrd":                   return (g, 1.2, g)
    if task in ("report_sample", "report_sample_TO"):
                                              return (g * 0.5, 1.2, g * 0.5)
    return (g, 1.2, g)

# Annotate every record
for rec in recs:
    stage, disp = classify(rec)
    rec["_stage"] = stage
    rec["_disp"]  = disp
    wm, rm, wallm = scale_for_rec(rec)
    rec["_work_mult"] = wm
    rec["_rss_mult"]  = rm
    rec["_wall_mult"] = wallm

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fmt_hms(s):
    if s is None: return "-"
    s = int(round(s))
    if s == 0: return "0s"
    d, r = divmod(s, 86400)
    h, r = divmod(r, 3600)
    m, sec = divmod(r, 60)
    if d:  return f"{d}d {h:02d}h {m:02d}m"
    if h:  return f"{h}h {m:02d}m {sec:02d}s"
    if m:  return f"{m}m {sec:02d}s"
    return f"{sec}s"

def fmt_gb(kb):
    if kb is None: return "-"
    gb = kb / 1024 / 1024
    return f"{gb/1024:.2f} TB" if gb >= 1024 else f"{gb:.2f} GB"

def fmt_pct(p):
    return "-" if p is None else f"{p:.0f}%"

def fmt_count(n):
    if n is None: return "-"
    if n >= 1e9: return f"{n/1e9:.2f}B"
    if n >= 1e6: return f"{n/1e6:.2f}M"
    if n >= 1e3: return f"{n/1e3:.1f}k"
    return str(int(n))

def h(s): return _html.escape(str(s))

# Sort key stored as data-val for JS sorter (raw numbers, so JS sorts numerically)
def td(val_fmt, sort_val, cls="num"):
    return f'<td class="{cls}" data-val="{sort_val}">{val_fmt}</td>'

def minibar(val, maxval, color="#4c78a8", width=90):
    pct = max(0.0, min(100.0, 100.0 * val / maxval)) if maxval else 0
    return (f'<div class="mbar"><div class="mbar-fill" '
            f'style="width:{pct:.1f}%;background:{color}"></div></div>')

# ---------------------------------------------------------------------------
# Aggregate structures
# ---------------------------------------------------------------------------
by_stage = defaultdict(list)
for rec in recs:
    by_stage[rec["_stage"]].append(rec)

# Scatter-collapse: group deepsomatic shards; one bucket per task display name
scatter_agg = {}
for rec in recs:
    key = "deepsomatic_paired (scatter)" if "deepsomatic_paired" in rec["call"] else rec["_disp"]
    if key not in scatter_agg:
        scatter_agg[key] = {
            "wall_s": 0, "shards": 0, "max_rss_kb": 0,
            "user_s_total": 0, "sys_s_total": 0, "cpu_pct_vals": [],
            "work_mult": 0, "rss_mult": 0, "wall_mult": 0,
            "stage": rec["_stage"], "disp": key,
        }
    b = scatter_agg[key]
    b["wall_s"]       += rec["wall_s"] or 0
    b["shards"]       += 1
    b["max_rss_kb"]    = max(b["max_rss_kb"], rec["max_rss_kb"] or 0)
    b["user_s_total"] += rec["user_s_total"] or 0
    b["sys_s_total"]  += rec["sys_s_total"] or 0
    if rec["cpu_pct_avg"] is not None:
        b["cpu_pct_vals"].append(rec["cpu_pct_avg"])
    b["work_mult"] = rec["_work_mult"]
    b["rss_mult"]  = rec["_rss_mult"]
    b["wall_mult"] = rec["_wall_mult"]

for b in scatter_agg.values():
    b["cpu_pct_avg"] = (sum(b["cpu_pct_vals"]) / len(b["cpu_pct_vals"])
                        if b["cpu_pct_vals"] else None)

# Stage-level totals (observed)
stage_obs = {}
for s, rs in by_stage.items():
    wall  = sum(r["wall_s"] or 0 for r in rs)
    user  = sum(r["user_s_total"] or 0 for r in rs)
    sys_  = sum(r["sys_s_total"]  or 0 for r in rs)
    peak  = max((r["max_rss_kb"] or 0 for r in rs), default=0)
    stage_obs[s] = {"wall": wall, "user": user, "sys": sys_, "peak": peak, "n": len(rs)}

# Stage-level WGS scaled
stage_wgs = {}
for s, rs in by_stage.items():
    wall = sum((r["wall_s"] or 0) * r["_work_mult"] for r in rs)
    user = sum((r["user_s_total"] or 0) * r["_work_mult"] for r in rs)
    sys_ = sum((r["sys_s_total"]  or 0) * r["_work_mult"] for r in rs)
    peak = max(((r["max_rss_kb"] or 0) * r["_rss_mult"] for r in rs), default=0)
    stage_wgs[s] = {"wall": wall, "user": user, "sys": sys_, "peak": peak}

# Overall totals
total_wall = sum(r["wall_s"] or 0 for r in recs)
total_user = sum(r["user_s_total"] or 0 for r in recs)
total_sys  = sum(r["sys_s_total"]  or 0 for r in recs)
total_cpu  = total_user + total_sys
peak_rss   = max((r["max_rss_kb"] or 0 for r in recs), default=0)

wgs_sum_wall = sum((r["wall_s"] or 0) * r["_work_mult"] for r in recs)
wgs_sum_user = sum((r["user_s_total"] or 0) * r["_work_mult"] for r in recs)
wgs_sum_sys  = sum((r["sys_s_total"]  or 0) * r["_work_mult"] for r in recs)
wgs_sum_cpu  = wgs_sum_user + wgs_sum_sys
wgs_peak_rss = max((r["max_rss_kb"] or 0) * r["_rss_mult"] for r in recs)

# WGS end-to-end wall estimate
def stage_crit_wall_wgs(s, rs):
    if s == "Small-variant calling":
        ds = [r for r in rs if "deepsomatic_paired" in r["call"]]
        other = [r for r in rs if "deepsomatic_paired" not in r["call"]]
        return max((r["wall_s"] or 0 for r in ds), default=0) * 3.0 + \
               sum((r["wall_s"] or 0) * r["_wall_mult"] for r in other)
    return sum((r["wall_s"] or 0) * r["_wall_mult"] for r in rs)

total_crit = sum(stage_crit_wall_wgs(s, by_stage[s]) for s in by_stage)
max_task_wgs = max((r["wall_s"] or 0) * r["_wall_mult"] for r in recs)
wgs_wall_est = max(total_crit / 1.8, max_task_wgs)

N_SAMPLES = 1000
thousand_cpu_hours  = wgs_sum_cpu  * N_SAMPLES / 3600
thousand_wall_hours = wgs_sum_wall * N_SAMPLES / 3600

BAM_PER_SAMPLE        = 200   # GB (tumor + normal ~30× HiFi)
storage_per_sample_gb = BAM_PER_SAMPLE * 3
storage_1000_tb       = storage_per_sample_gb * N_SAMPLES / 1024

# Workflow timing from log
log_lines = open(os.path.join(ROOT, "workflow.log")).readlines()
def _ts(line):
    m = _re.match(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)", line)
    return datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S.%f") if m else None
wf_start = _ts(log_lines[0])
wf_end   = None
for ln in reversed(log_lines):
    t = _ts(ln)
    if t: wf_end = t; break
clock_elapsed = (wf_end - wf_start).total_seconds() if (wf_start and wf_end) else None

# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------

def sortable_table(table_id, col_headers, rows_html, caption=None):
    """Wrap rows in a full <table> with sort-enabled <th> elements."""
    ths = "".join(
        f'<th class="sortable" onclick="sortTable(\'{table_id}\',{i})">'
        f'{h(c)} <span class="sort-icon">⇅</span></th>'
        for i, c in enumerate(col_headers)
    )
    cap = f"<caption>{h(caption)}</caption>" if caption else ""
    return f"""<table class="data sortable-tbl" id="{table_id}">
  {cap}
  <thead><tr>{ths}</tr></thead>
  <tbody>{rows_html}</tbody>
</table>"""

# ── Stage summary table ──────────────────────────────────────────────────────
def stage_summary_rows():
    max_user_obs = max((stage_obs[s]["user"] for s in stage_obs), default=1) or 1
    max_wall_obs = max((stage_obs[s]["wall"] for s in stage_obs), default=1) or 1
    max_rss_obs  = max((stage_obs[s]["peak"] for s in stage_obs), default=1) or 1
    max_user_wgs = max((stage_wgs[s]["user"] for s in stage_wgs), default=1) or 1
    rows = []
    for s in STAGE_ORDER:
        if s not in stage_obs: continue
        o = stage_obs[s]; w = stage_wgs[s]
        color = STAGE_COLOR.get(s, "#aaa")
        cpu_pct_obs = (o["user"] + o["sys"]) / o["wall"] * 100 if o["wall"] else 0
        rows.append(
            f"<tr>"
            f'<td><span class="dot" style="background:{color}"></span>{h(s)}</td>'
            + td(str(o["n"]),            o["n"])
            + td(fmt_hms(o["wall"]),     o["wall"])
            + td(fmt_hms(o["user"]),     o["user"],  "num")
            + f'<td class="num bar-cell" data-val="{o["user"]}">'
              f'{minibar(o["user"], max_user_obs, color)}</td>'
            + td(fmt_pct(cpu_pct_obs),   cpu_pct_obs)
            + td(fmt_gb(o["peak"]),      o["peak"])
            + f'<td class="num bar-cell" data-val="{o["peak"]}">'
              f'{minibar(o["peak"], max_rss_obs, color)}</td>'
            + td(fmt_hms(w["wall"]),     w["wall"])
            + td(fmt_hms(w["user"]),     w["user"])
            + f'<td class="num bar-cell" data-val="{w["user"]}">'
              f'{minibar(w["user"], max_user_wgs, color)}</td>'
            + td(fmt_hms(w["user"] * N_SAMPLES), w["user"] * N_SAMPLES)
            + f"</tr>\n"
        )
    return "".join(rows)

# ── Task-level table (scatter-collapsed) ─────────────────────────────────────
def task_rows_all():
    tasks = sorted(scatter_agg.values(), key=lambda b: b["user_s_total"], reverse=True)
    max_user = max((b["user_s_total"] for b in tasks), default=1) or 1
    max_wall = max((b["wall_s"]       for b in tasks), default=1) or 1
    max_rss  = max((b["max_rss_kb"]   for b in tasks), default=1) or 1
    rows = []
    for b in tasks:
        color = STAGE_COLOR.get(b["stage"], "#aaa")
        wgs_user = b["user_s_total"] * b["work_mult"]
        wgs_peak = b["max_rss_kb"]   * b["rss_mult"]
        cpu_eff  = (b["user_s_total"] + b["sys_s_total"]) / b["wall_s"] * 100 \
                   if b["wall_s"] else 0
        rows.append(
            f"<tr>"
            f'<td><span class="dot" style="background:{color}"></span>{h(b["stage"])}</td>'
            f'<td>{h(b["disp"])}</td>'
            + td(str(b["shards"]), b["shards"])
            + td(fmt_hms(b["wall_s"]),      b["wall_s"])
            + f'<td class="num bar-cell" data-val="{b["wall_s"]}">'
              f'{minibar(b["wall_s"],       max_wall, color)}</td>'
            + td(fmt_hms(b["user_s_total"]), b["user_s_total"])
            + f'<td class="num bar-cell" data-val="{b["user_s_total"]}">'
              f'{minibar(b["user_s_total"], max_user, color)}</td>'
            + td(fmt_pct(cpu_eff),           cpu_eff)
            + td(fmt_gb(b["max_rss_kb"]),    b["max_rss_kb"])
            + f'<td class="num bar-cell" data-val="{b["max_rss_kb"]}">'
              f'{minibar(b["max_rss_kb"],    max_rss, color)}</td>'
            + td(fmt_hms(wgs_user),          wgs_user)
            + td(fmt_gb(wgs_peak),           wgs_peak)
            + f"</tr>\n"
        )
    return "".join(rows)

# ── Per-stage drilldown ───────────────────────────────────────────────────────
def per_stage_drilldown():
    blocks = []
    for s in STAGE_ORDER:
        if s not in by_stage: continue
        rs = sorted(by_stage[s], key=lambda r: r["user_s_total"] or 0, reverse=True)
        color = STAGE_COLOR.get(s, "#aaa")

        # Separate deepsomatic shards from others
        ds_shards = [r for r in rs if "deepsomatic_paired" in r["call"]]
        other     = [r for r in rs if "deepsomatic_paired" not in r["call"]]

        max_user_s = max((r["user_s_total"] or 0 for r in rs), default=1) or 1
        max_rss_s  = max((r["max_rss_kb"]   or 0 for r in rs), default=1) or 1

        tbl_rows = []

        # Individual (non-scatter) task rows
        for r in other:
            cpu_eff = ((r["user_s_total"] or 0) + (r["sys_s_total"] or 0)) / r["wall_s"] * 100 \
                      if r["wall_s"] else 0
            tbl_rows.append(
                f"<tr>"
                f"<td><code>{h(r['_disp'])}</code></td>"
                + td(fmt_hms(r["wall_s"]),      r["wall_s"])
                + td(fmt_hms(r["user_s_total"]), r["user_s_total"] or 0)
                + f'<td class="num bar-cell" data-val="{r["user_s_total"] or 0}">'
                  f'{minibar(r["user_s_total"] or 0, max_user_s, color)}</td>'
                + td(fmt_pct(cpu_eff),           cpu_eff)
                + td(fmt_gb(r["max_rss_kb"]),    r["max_rss_kb"] or 0)
                + f'<td class="num bar-cell" data-val="{r["max_rss_kb"] or 0}">'
                  f'{minibar(r["max_rss_kb"] or 0, max_rss_s, color)}</td>'
                + td(fmt_hms(r["sys_s_total"]),  r["sys_s_total"] or 0)
                + td(str(r["n_iters"]),           r["n_iters"])
                + f"</tr>\n"
            )

        # Scatter summary rows
        if ds_shards:
            n         = len(ds_shards)
            tot_wall  = sum(r["wall_s"] or 0 for r in ds_shards)
            max_wall  = max(r["wall_s"] or 0 for r in ds_shards)
            min_wall  = min(r["wall_s"] or 0 for r in ds_shards)
            tot_user  = sum(r["user_s_total"] or 0 for r in ds_shards)
            tot_sys   = sum(r["sys_s_total"]  or 0 for r in ds_shards)
            pk_rss    = max(r["max_rss_kb"]   or 0 for r in ds_shards)
            mean_cpu  = sum(r["cpu_pct_avg"]  or 0 for r in ds_shards) / n
            cpu_eff   = (tot_user + tot_sys) / tot_wall * 100 if tot_wall else 0
            tbl_rows.append(
                f"<tr class='scatter-sum'>"
                f"<td><code>deepsomatic_paired ({n} shards — sum)</code></td>"
                + td(fmt_hms(tot_wall), tot_wall)
                + td(fmt_hms(tot_user), tot_user)
                + f'<td class="num bar-cell" data-val="{tot_user}">'
                  f'{minibar(tot_user, max_user_s, color)}</td>'
                + td(fmt_pct(cpu_eff), cpu_eff)
                + td(fmt_gb(pk_rss),   pk_rss)
                + f'<td class="num bar-cell" data-val="{pk_rss}">'
                  f'{minibar(pk_rss, max_rss_s, color)}</td>'
                + td(fmt_hms(tot_sys), tot_sys)
                + td(str(n), n)
                + f"</tr>\n"
            )
            tbl_rows.append(
                f"<tr class='scatter-sum'>"
                f"<td><code>deepsomatic_paired per-shard wall: "
                f"min {fmt_hms(min_wall)} / max {fmt_hms(max_wall)}</code></td>"
                f'<td colspan="8" class="num muted">'
                f"avg CPU {fmt_pct(mean_cpu)}, peak RSS {fmt_gb(pk_rss)}</td>"
                f"</tr>\n"
            )

        tid = f"stg_{s.replace(' ','_').replace('/','_').replace('-','_')}"
        tbl = sortable_table(
            tid,
            ["Task", "Wall", "User CPU", "CPU bar", "CPU eff%", "Peak RSS", "RSS bar", "Sys CPU", "Iter"],
            "".join(tbl_rows),
        )

        o = stage_obs[s]
        w = stage_wgs[s]
        kpis = (
            f'<div class="stage-kpis">'
            f'<span class="skpi"><b>{o["n"]}</b> calls</span>'
            f'<span class="skpi"><b>{fmt_hms(o["wall"])}</b> sum-wall</span>'
            f'<span class="skpi"><b>{fmt_hms(o["user"])}</b> user CPU</span>'
            f'<span class="skpi"><b>{fmt_gb(o["peak"])}</b> peak RSS</span>'
            f'<span class="skpi muted">WGS est: <b>{fmt_hms(w["user"])}</b> user CPU</span>'
            f'</div>'
        )
        blocks.append(
            f'<details class="stage-detail" id="detail_{s.replace(" ","_")}">\n'
            f'  <summary><span class="dot lg" style="background:{color}"></span>'
            f'  <strong>{h(s)}</strong>{kpis}</summary>\n'
            f'  {tbl}\n'
            f'</details>\n'
        )
    return "\n".join(blocks)

# ── Individual task call table (full, unsummarised) ─────────────────────────
def all_task_calls_rows():
    rs = sorted(recs, key=lambda r: r["user_s_total"] or 0, reverse=True)
    max_user = max((r["user_s_total"] or 0 for r in rs), default=1) or 1
    rows = []
    for r in rs:
        color = STAGE_COLOR.get(r["_stage"], "#aaa")
        cpu_eff = ((r["user_s_total"] or 0) + (r["sys_s_total"] or 0)) / r["wall_s"] * 100 \
                  if r["wall_s"] else 0
        rows.append(
            f"<tr>"
            f'<td><span class="dot" style="background:{color}"></span>{h(r["_stage"])}</td>'
            f'<td>{h(r["_disp"])}</td>'
            f'<td><code class="smallcode">{h(r["call"])}</code></td>'
            + td(fmt_hms(r["wall_s"]),       r["wall_s"])
            + td(fmt_hms(r["user_s_total"]), r["user_s_total"] or 0)
            + f'<td class="num bar-cell" data-val="{r["user_s_total"] or 0}">'
              f'{minibar(r["user_s_total"] or 0, max_user, color)}</td>'
            + td(fmt_pct(cpu_eff),           cpu_eff)
            + td(fmt_gb(r["max_rss_kb"]),    r["max_rss_kb"] or 0)
            + td(fmt_hms(r["sys_s_total"]),  r["sys_s_total"] or 0)
            + f"</tr>\n"
        )
    return "".join(rows)

# ── Observations ────────────────────────────────────────────────────────────
rss_ranked  = sorted((r for r in recs if r["max_rss_kb"]),
                     key=lambda r: r["max_rss_kb"], reverse=True)[:5]
slow_single = sorted(
    [r for r in recs if r["wall_s"] and r["wall_s"] > 60
     and r["cpu_pct_avg"] is not None and r["cpu_pct_avg"] < 120],
    key=lambda r: r["wall_s"], reverse=True
)[:6]

ds_recs         = [r for r in recs if "deepsomatic_paired" in r["call"]]
phase_n_recs    = [r for r in recs if "phaseNormalBam"    in r["call"]]
align_recs      = [r for r in recs if "align_single_bam"  in r["call"]]
ds_max_wgs      = max((r["wall_s"] or 0) * r["_wall_mult"] for r in ds_recs)      if ds_recs      else 0
phase_n_max_wgs = max((r["wall_s"] or 0) * r["_wall_mult"] for r in phase_n_recs) if phase_n_recs else 0
align_max_wgs   = max((r["wall_s"] or 0) * r["_wall_mult"] for r in align_recs)   if align_recs   else 0

# ============================================================================
# HTML
# ============================================================================
html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>HiFi-somatic-WDL Resource Utilization Report</title>
<style>
:root {{
  --fg:#1f2933; --muted:#52606d; --accent:#4c78a8; --accent2:#f58518;
  --rule:#cbd2d9; --bg:#f5f7f9; --card:#ffffff; --thead:#e4e7eb; --scatter:#fff8e7;
  --radius:6px;
}}
* {{ box-sizing:border-box; margin:0; padding:0; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
        color:var(--fg); background:var(--bg); padding:1.5rem 2rem; line-height:1.5; }}
.container {{ max-width:1320px; margin:0 auto; }}
h1 {{ font-size:1.8rem; margin-bottom:.2rem; }}
h2 {{ font-size:1.2rem; margin:2rem 0 .6rem;
      border-bottom:3px solid var(--accent); padding-bottom:.3rem; }}
h3 {{ font-size:1rem; margin:1.2rem 0 .4rem; color:var(--muted); }}
.meta {{ color:var(--muted); font-size:.88rem; margin-bottom:1.2rem; }}
.meta b {{ color:var(--fg); }}

/* KPI tiles */
.kpi-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
             gap:.6rem; margin:.8rem 0 1.4rem; }}
.kpi {{ background:var(--card); border:1px solid var(--rule); border-radius:var(--radius);
        padding:.65rem 1rem; }}
.kpi .val {{ font-size:1.2rem; font-weight:700; color:var(--accent); }}
.kpi .val.warn {{ color:var(--accent2); }}
.kpi .lbl {{ color:var(--muted); font-size:.8rem; margin-top:.1rem; }}

/* Tables */
table.data {{ border-collapse:collapse; width:100%; background:var(--card);
              font-size:.83rem; margin:.4rem 0 1rem; border-radius:var(--radius);
              overflow:hidden; box-shadow:0 1px 4px rgba(0,0,0,.06); }}
table.data th, table.data td {{ padding:.42rem .65rem; border-bottom:1px solid var(--rule);
                                  text-align:left; vertical-align:middle; }}
table.data thead {{ background:var(--thead); }}
th.sortable {{ cursor:pointer; user-select:none; white-space:nowrap; }}
th.sortable:hover {{ background:#d5dce5; }}
.sort-icon {{ opacity:.45; font-size:.8em; }}
th.sort-asc .sort-icon::after  {{ content:"▲"; opacity:1; }}
th.sort-desc .sort-icon::after {{ content:"▼"; opacity:1; }}
th.sort-asc  .sort-icon, th.sort-desc .sort-icon {{ opacity:1; }}
td.num {{ font-variant-numeric:tabular-nums; white-space:nowrap; }}
td.bar-cell {{ min-width:100px; }}
tr.scatter-sum {{ background:var(--scatter); }}
tr:hover td {{ background:#f0f4f8; }}

/* Mini inline bar */
.mbar {{ display:inline-block; width:88px; height:.55rem; background:#e9ecef;
         border-radius:3px; vertical-align:middle; }}
.mbar-fill {{ height:100%; border-radius:3px; transition:width .2s; }}

code {{ font-family:SFMono-Regular,Consolas,Menlo,monospace; font-size:.78rem;
        background:#eef1f3; padding:.1rem .3rem; border-radius:3px; white-space:nowrap; }}
.smallcode {{ font-size:.72rem; }}

/* Stage dots */
.dot {{ display:inline-block; width:.65rem; height:.65rem; border-radius:50%;
        margin-right:.35rem; vertical-align:middle; flex-shrink:0; }}
.dot.lg {{ width:.9rem; height:.9rem; margin-right:.5rem; }}

/* Stage detail accordion */
details.stage-detail {{ background:var(--card); border:1px solid var(--rule);
                         border-radius:var(--radius); margin:.5rem 0;
                         box-shadow:0 1px 3px rgba(0,0,0,.04); }}
details.stage-detail > summary {{ list-style:none; padding:.7rem 1rem;
    display:flex; align-items:center; cursor:pointer; }}
details.stage-detail > summary::-webkit-details-marker {{ display:none; }}
details.stage-detail[open] > summary {{ border-bottom:1px solid var(--rule); }}
details.stage-detail > summary:hover {{ background:#f7f9fb; }}
.stage-kpis {{ display:flex; gap:.9rem; margin-left:auto; flex-wrap:wrap; }}
.skpi {{ font-size:.82rem; color:var(--muted); }}
.skpi b {{ color:var(--fg); }}
.skpi.muted {{ opacity:.7; }}
details.stage-detail table {{ margin:.6rem .8rem .8rem; width:calc(100% - 1.6rem); }}

/* Callouts */
.callout {{ border-left:4px solid var(--accent); background:#f3f7fb;
            padding:.7rem 1rem; margin:.7rem 0; border-radius:0 var(--radius) var(--radius) 0; }}
.callout.warn {{ border-left-color:var(--accent2); background:#fff7ec; }}
ul.obs {{ padding-left:1.2rem; margin:.4rem 0; }}
ul.obs li {{ margin-bottom:.25rem; font-size:.88rem; }}
.muted {{ color:var(--muted); }}

footer {{ color:var(--muted); font-size:.82rem; margin-top:2rem;
          padding-top:.8rem; border-top:1px solid var(--rule); }}

/* TOC nav */
nav.toc {{ background:var(--card); border:1px solid var(--rule); border-radius:var(--radius);
           padding:.6rem 1rem; margin-bottom:1.5rem; font-size:.85rem; }}
nav.toc a {{ color:var(--accent); text-decoration:none; margin-right:1rem; }}
nav.toc a:hover {{ text-decoration:underline; }}
</style>
</head>
<body>
<div class="container">

<h1>HiFi-somatic-WDL &mdash; Resource Utilization Report</h1>
<div class="meta">
  <b>Run:</b> <code>{h(ROOT)}</code> &nbsp;|&nbsp;
  <b>Start:</b> {h(wf_start.strftime('%Y-%m-%d %H:%M:%S')) if wf_start else 'n/a'} &nbsp;|&nbsp;
  <b>End:</b>   {h(wf_end.strftime('%Y-%m-%d %H:%M:%S'))   if wf_end   else 'n/a'} &nbsp;|&nbsp;
  <b>Elapsed:</b> {h(fmt_hms(clock_elapsed))} &nbsp;|&nbsp;
  <b>Tasks profiled:</b> {len(recs)} &nbsp;|&nbsp;
  <b>Generated:</b> {h(datetime.now().strftime('%Y-%m-%d %H:%M'))}
</div>

<nav class="toc">
  <b>Jump to:</b>
  <a href="#sec-headline">Headline metrics</a>
  <a href="#sec-stages">Resource by stage</a>
  <a href="#sec-tasks">Resource by task</a>
  <a href="#sec-drilldown">Stage drilldown</a>
  <a href="#sec-allcalls">All task calls</a>
  <a href="#sec-wgs">WGS projections</a>
  <a href="#sec-obs">Observations</a>
  <a href="#sec-method">Methodology</a>
</nav>

<!-- ═══════════════════════════════════════════════════════ §1 HEADLINE ══ -->
<h2 id="sec-headline">1 &mdash; Headline metrics (chr22 observed)</h2>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(recs)}</div>
    <div class="lbl">task calls profiled</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(clock_elapsed))}</div>
    <div class="lbl">workflow wall-clock</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(total_wall))}</div>
    <div class="lbl">sum of per-task wall</div></div>
  <div class="kpi"><div class="val">{total_wall/clock_elapsed:.1f}×</div>
    <div class="lbl">effective parallelism</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(total_user))}</div>
    <div class="lbl">total user CPU</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(total_sys))}</div>
    <div class="lbl">total system CPU</div></div>
  <div class="kpi"><div class="val warn">{h(fmt_gb(peak_rss))}</div>
    <div class="lbl">peak max-RSS (any task)</div></div>
  <div class="kpi"><div class="val">{(total_cpu/total_wall*100) if total_wall else 0:.0f}%</div>
    <div class="lbl">CPU efficiency (user+sys ÷ sum-wall)</div></div>
</div>

<!-- ═══════════════════════════════════════════════════════ §2 STAGES ════ -->
<h2 id="sec-stages">2 &mdash; Resource utilization by stage</h2>
<p class="muted" style="font-size:.85rem;margin-bottom:.5rem">
  Click any column header to sort. Bar columns reflect relative magnitude within
  the observed chr22 run. <em>WGS user CPU</em> and <em>×1000 CPU</em> are scaled
  projections — see §7 for methodology.
</p>
{sortable_table(
    "tbl_stages",
    ["Stage","Calls","Sum wall","User CPU","CPU bar","CPU eff%","Peak RSS","RSS bar",
     "WGS sum-wall","WGS user CPU","WGS CPU bar","×1000 user CPU"],
    stage_summary_rows()
)}

<!-- ═══════════════════════════════════════════════════════ §3 TASKS ═════ -->
<h2 id="sec-tasks">3 &mdash; Resource utilization by task (scatter collapsed)</h2>
<p class="muted" style="font-size:.85rem;margin-bottom:.5rem">
  Deepsomatic scatter shards are aggregated into one row. Sorted by user CPU descending
  by default. Click any header to re-sort.
</p>
{sortable_table(
    "tbl_tasks",
    ["Stage","Task","Shards","Sum wall","Wall bar","User CPU","CPU bar","CPU eff%",
     "Peak RSS","RSS bar","WGS user CPU","WGS peak RSS"],
    task_rows_all()
)}

<!-- ═══════════════════════════════════════════════════════ §4 DRILLDOWN ═ -->
<h2 id="sec-drilldown">4 &mdash; Stage drilldown</h2>
<p class="muted" style="font-size:.85rem;margin-bottom:.6rem">
  Click a stage to expand its task-level detail. Tables within each stage are also sortable.
</p>
{per_stage_drilldown()}

<!-- ═══════════════════════════════════════════════════════ §5 ALL CALLS ═ -->
<h2 id="sec-allcalls">5 &mdash; All individual task calls</h2>
<p class="muted" style="font-size:.85rem;margin-bottom:.5rem">
  Every profiled call, unsummarised. Sorted by user CPU descending.
</p>
{sortable_table(
    "tbl_all",
    ["Stage","Task","Call path","Wall","User CPU","CPU bar","CPU eff%","Peak RSS","Sys CPU"],
    all_task_calls_rows()
)}

<!-- ═══════════════════════════════════════════════════════ §6 WGS ════════ -->
<h2 id="sec-wgs">6 &mdash; WGS scaling projections</h2>

<div class="callout">
  <strong>Scaling basis.</strong>
  chr22 = {fmt_count(CHR22_BP)} bp; GRCh38 primary assembly = {fmt_count(WGS_BP)} bp →
  genome-to-chr22 ratio ≈ <strong>{GENOME_SCALE:.1f}×</strong>.
  CPU and wall scale linearly for data-bound tools; peak RSS assumed ×1.0–1.5
  (reference-index dominated). <code>deepsomatic_paired</code> per-shard wall uses 3×
  (chr1 ≈ 3× chr22 bp).
</div>

<h3>Single WGS sample</h3>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{h(fmt_hms(wgs_sum_wall))}</div>
    <div class="lbl">projected sum-of-task wall</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(wgs_wall_est))}</div>
    <div class="lbl">projected end-to-end wall (÷1.8× parallelism)</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(wgs_sum_user))}</div>
    <div class="lbl">projected user CPU</div></div>
  <div class="kpi"><div class="val warn">{h(fmt_gb(wgs_peak_rss))}</div>
    <div class="lbl">projected peak RSS</div></div>
  <div class="kpi"><div class="val">{BAM_PER_SAMPLE} GB</div>
    <div class="lbl">input BAM footprint (T+N, 30× HiFi)</div></div>
  <div class="kpi"><div class="val">≈ {storage_per_sample_gb} GB</div>
    <div class="lbl">intermediate storage est.</div></div>
</div>

<h3>1 000 WGS samples (same cluster)</h3>
<div class="kpi-grid">
  <div class="kpi"><div class="val warn">{thousand_cpu_hours:,.0f}</div>
    <div class="lbl">total CPU-hours</div></div>
  <div class="kpi"><div class="val">{thousand_wall_hours:,.0f} hr</div>
    <div class="lbl">sum-of-task wall (serialised)</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(wgs_wall_est * N_SAMPLES / 100))}</div>
    <div class="lbl">wall @ 100 samples concurrent</div></div>
  <div class="kpi"><div class="val">{h(fmt_hms(wgs_wall_est))}</div>
    <div class="lbl">wall @ 1 000 samples concurrent</div></div>
  <div class="kpi"><div class="val">{storage_1000_tb:,.0f} TB</div>
    <div class="lbl">total storage estimate</div></div>
</div>

<!-- ═══════════════════════════════════════════════════════ §7 OBS ════════ -->
<h2 id="sec-obs">7 &mdash; Observations</h2>

<div class="callout">
  <strong>Top memory consumers (chr22 observed):</strong>
  <ul class="obs">
    {''.join(f"<li><code>{h(r['call'])}</code> &mdash; {h(fmt_gb(r['max_rss_kb']))} ({h(r['_disp'])})</li>" for r in rss_ranked)}
  </ul>
</div>

<div class="callout warn">
  <strong>Long-running low-CPU tasks (wall &gt; 1 min, CPU% &lt; 120) — over-provisioned cores:</strong>
  <ul class="obs">
    {''.join(f"<li><code>{h(r['call'])}</code> &mdash; wall {h(fmt_hms(r['wall_s']))}, CPU {h(fmt_pct(r['cpu_pct_avg']))} ({h(r['_disp'])})</li>" for r in slow_single)}
  </ul>
  These hold multiple cores in the <code>runtime</code> block but run single-threaded.
  At WGS scale wall grows ≈{GENOME_SCALE:.0f}×, magnifying the idle-core waste.
</div>

<div class="callout">
  <strong>Critical-path bottlenecks at WGS scale:</strong>
  <ul class="obs">
    <li><code>deepsomatic_paired</code> heaviest-shard projected wall ≈ <b>{h(fmt_hms(ds_max_wgs))}</b>
        — drives small-variant stage critical path.</li>
    <li><code>phaseNormalBam</code> (hiphase) projected wall ≈ <b>{h(fmt_hms(phase_n_max_wgs))}</b>
        — single-task stage, its wall <em>is</em> the phasing critical path.</li>
    <li><code>align_single_bam</code> (pbmm2) projected wall ≈ <b>{h(fmt_hms(align_max_wgs))}</b>
        — also the memory floor; each pbmm2 instance requires &gt;90 GB RAM.</li>
  </ul>
</div>

<!-- ═══════════════════════════════════════════════════════ §8 METHOD ════ -->
<h2 id="sec-method">8 &mdash; Methodology</h2>
<ul class="obs">
  <li>Metrics sourced from each task's <code>out/resource_usage/resource_usage.tsv</code>,
      produced by <code>/usr/bin/time -v</code> wrappers in the WDL tasks.</li>
  <li><code>wall_clock_seconds</code> covers the entire <code>command&nbsp;{{}}&nbsp;{{}}</code>
      block including pre/post-processing; the timed command wall may be slightly shorter.</li>
  <li>For <code>prioritize_dmr_intogen</code> (multi-iteration append mode), user/sys CPU
      are summed across iterations; max RSS and mean CPU% are taken.</li>
  <li><strong>Scatter collapse:</strong> <code>deepsomatic_paired</code> shards are
      merged into one row (sum wall/CPU, max RSS) for the task-level table.</li>
  <li><strong>CPU efficiency</strong> = (user + sys) ÷ sum-wall × 100; 100% means
      every allocated core was busy at all times.</li>
  <li><strong>Effective parallelism</strong> = sum-of-task-wall ÷ workflow-elapsed;
      values &gt;1 reflect concurrent task execution.</li>
  <li><strong>WGS scaling:</strong> genome-size ratio chr22→GRCh38 = {GENOME_SCALE:.2f}×;
      data-bound tasks scale linearly; RSS ×1.0–1.5 (index-dominated);
      <code>deepsomatic_paired</code> per-shard wall ×3 (chr1 ≈ 3× chr22).</li>
</ul>

<footer>
  Generated {h(datetime.now().strftime('%Y-%m-%d %H:%M'))} from
  <code>{h(ROOT)}</code> &mdash; HiFi-somatic-WDL resource profiling.
</footer>
</div>

<script>
// ── Table sorter ──────────────────────────────────────────────────────────
function sortTable(tableId, colIdx) {{
  const tbl  = document.getElementById(tableId);
  const ths  = tbl.querySelectorAll('thead th');
  const th   = ths[colIdx];
  const tbody = tbl.querySelector('tbody');
  const rows  = Array.from(tbody.querySelectorAll('tr'));

  // Determine sort direction
  let asc = true;
  if (th.classList.contains('sort-asc'))  {{ asc = false; }}
  if (th.classList.contains('sort-desc')) {{ asc = true;  }}

  // Clear all sort indicators
  ths.forEach(t => t.classList.remove('sort-asc','sort-desc'));
  th.classList.add(asc ? 'sort-asc' : 'sort-desc');

  rows.sort((a, b) => {{
    const tds_a = a.querySelectorAll('td');
    const tds_b = b.querySelectorAll('td');
    if (colIdx >= tds_a.length || colIdx >= tds_b.length) return 0;

    // Use data-val if present (raw number for proper numeric sort)
    const raw_a = tds_a[colIdx].dataset.val;
    const raw_b = tds_b[colIdx].dataset.val;

    let va, vb;
    if (raw_a !== undefined && raw_b !== undefined) {{
      va = parseFloat(raw_a);
      vb = parseFloat(raw_b);
      if (isNaN(va)) va = -Infinity;
      if (isNaN(vb)) vb = -Infinity;
    }} else {{
      va = tds_a[colIdx].textContent.trim().toLowerCase();
      vb = tds_b[colIdx].textContent.trim().toLowerCase();
      // Try numeric parse for text cells too
      const nva = parseFloat(va), nvb = parseFloat(vb);
      if (!isNaN(nva) && !isNaN(nvb)) {{ va = nva; vb = nvb; }}
    }}

    if (va < vb) return asc ?  -1 : 1;
    if (va > vb) return asc ?   1 : -1;
    return 0;
  }});

  rows.forEach(r => tbody.appendChild(r));
}}

// Open the stage detail section whose anchor is in the URL hash
(function() {{
  const hash = window.location.hash;
  if (hash && hash.startsWith('#detail_')) {{
    const el = document.querySelector(hash);
    if (el) el.open = true;
  }}
}})();
</script>
</body>
</html>
"""

with open(OUT_HTML, "w") as f:
    f.write(html)

print("Wrote:", OUT_HTML)
print(f"\nObserved chr22 run:")
print(f"  elapsed              = {fmt_hms(clock_elapsed)}")
print(f"  sum-of-task wall     = {fmt_hms(total_wall)}")
print(f"  effective parallelism= {total_wall/clock_elapsed:.1f}x")
print(f"  total user CPU       = {fmt_hms(total_user)}")
print(f"  total sys CPU        = {fmt_hms(total_sys)}")
print(f"  peak RSS             = {fmt_gb(peak_rss)}")
print(f"\nSingle WGS projections:")
print(f"  sum-wall             = {fmt_hms(wgs_sum_wall)}")
print(f"  end-to-end wall est  = {fmt_hms(wgs_wall_est)}")
print(f"  total CPU            = {fmt_hms(wgs_sum_cpu)}")
print(f"  peak RSS             = {fmt_gb(wgs_peak_rss)}")
print(f"\n1000 WGS:")
print(f"  CPU-hours            = {thousand_cpu_hours:,.0f}")
print(f"  wall-hours serial    = {thousand_wall_hours:,.0f}")
