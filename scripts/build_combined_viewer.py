"""
Build the unified SR manifold viewer (v2).
Single right-side filter panel with four collapsible sections:
  1. Manifold  — dropdown: All / Surface / Ridge / ω / Surface+Ridge / None
  2. Other species (Naveh) — per-species checkboxes with emoji
  3. NIA ITP — intervention dropdown + sex
  4. ZIMS zoo animals — class / sex / arm / RMS / genus
Three tabs: 3D Manifold | 2D Analysis | Phylo Analysis
"""
import pandas as pd, numpy as np, json, math

ITP_CSV  = (
    "/Users/Avimayo/Library/CloudStorage/GoogleDrive-mayoavi@gmail.com"
    "/.shortcut-targets-by-id/1JUcTeGwLojxa9QVWv-SdPTq9kAz4r-xR"
    "/SysMedBook2021/Chapter 6_Aging/ITP survival analysis/itp_fp_fits.csv"
)
ZIMS_CSV = "/Users/Avimayo/Documents/saturated-removal-manifold/results/zims_all.csv"
TAX_CSV  = "/Users/Avimayo/Documents/saturated-removal-manifold/results/zims_taxonomy.csv"
HTML_IN  = "/Users/Avimayo/Documents/saturated-removal-manifold/index.html"
HTML_OUT = "/Users/Avimayo/Documents/saturated-removal-manifold/viewer.html"
DIV_ID   = "a1d89836-0daf-4582-ae54-c7c4aca77de6"

# ── Manifold trace map (from index.html inspection) ───────────────────────────
MANIFOLD_TRACES = [
    (0, "constraint surface", "Surface",            "🗺"),
    (1, "model ridge (Eq.)",  "Ridge",              "〰"),
    (2, "ω on ridge",         "ω labels",           "Ω"),
]
NAVEH_SPECIES = [
    (3,  "Mice",                    "🐭"),
    (4,  "Yeast",                   "🧫"),
    (5,  "E. coli",                 "🦠"),
    (6,  "Cats",                    "🐱"),
    (7,  "Drosophila",              "🪰"),
    (8,  "Dogs",                    "🐶"),
    (9,  "C. elegans",              "🪱"),
    (10, "Guinea pig",              "🐹"),
    (11, "Humans",                  "👤"),
    (12, "Killifish",               "🐟"),
    (13, "individual measurements", "📊"),
]
NAVEH_IDX = [s[0] for s in NAVEH_SPECIES]

# ══════════════════════════════════════════════════════════════════════════════
# ITP data
# ══════════════════════════════════════════════════════════════════════════════
itp = pd.read_csv(ITP_CSV)
itp["lx"] = np.log10(itp["rho_beta"])
itp["ly"] = np.log10(itp["rho_eta"])
itp["lz"] = np.log10(itp["rho_eps"])

SEX_COLOR  = {"f": "#e8634a", "m": "#2c7bb6"}
ARM_SYMBOL = {"removal": "circle", "production": "diamond"}

ITP_CATEGORIES = {
    "Rapamycin":    ["Rapa_hi","Rapa_lo","Rapa_mid",
                     "Rapa_hi_continuous","Rapa_hi_cycle","Rapa_hi_start_stop"],
    "17α-E2":       ["17aE2","17aE2_hi","17aE2_16m","17aE2_20m"],
    "Metabolic":    ["Met","MetRapa","ACA","Cana"],
    "Antioxidants": ["MitoQ","GGA","MB","bGPA","NR","Min","CC"],
    "Other":        ["HBX","INT-767","DMAG","MIF098","Prot","UDCA"],
    "Control":      ["Control"],
}

itp_traces = []
for grp in sorted(itp["group"].unique()):
    for sex in ["f","m"]:
        sub = itp[(itp["group"]==grp)&(itp["sex"]==sex)]
        if sub.empty: continue
        is_ctrl = grp == "Control"
        color   = SEX_COLOR[sex]
        hover = [
            f"<b>{r['group']}</b> | {r['cohort']} | {'♀' if r['sex']=='f' else '♂'}<br>"
            f"L={r['L']:.0f} d  s={r['s']:.2f}  rms={r['rms']:.3f}  κ={r['kappa']:.4f}<br>"
            f"arm={r['arm']}<br>"
            f"ρ<sub>β</sub>={r['rho_beta']:.3f}  ρ<sub>η</sub>={r['rho_eta']:.3f}  ρ<sub>ε</sub>={r['rho_eps']:.3f}"
            for _,r in sub.iterrows()
        ]
        itp_traces.append({
            "type":"scatter3d",
            "name":f"{grp} ({'♀' if sex=='f' else '♂'})",
            "x":sub["lx"].tolist(),"y":sub["ly"].tolist(),"z":sub["lz"].tolist(),
            "text":hover,"hoverinfo":"text","mode":"markers",
            "marker":{
                "color":color,
                "symbol":[ARM_SYMBOL[a] for a in sub["arm"]],
                "size":7 if is_ctrl else 5,"opacity":0.95,
                "line":{"color":"black" if is_ctrl else color,
                        "width":1.5 if is_ctrl else 0},
            },
            "legendgroup":grp,"showlegend":True,
            "meta":{"group":grp,"sex":sex,"sex_color":color,
                    "kappa":sub["kappa"].tolist(),"dataset":"itp"},
        })

itp_opts  = '<option value="all">— All interventions —</option>\n'
for cat, members in ITP_CATEGORIES.items():
    itp_opts += f'<optgroup label="── {cat} ──">\n'
    itp_opts += f'  <option value="cat:{cat}">★ All {cat}</option>\n'
    for g in members:
        itp_opts += f'  <option value="{g}">&nbsp;&nbsp;{g}</option>\n'
    itp_opts += '</optgroup>\n'
itp_cat_map_json = json.dumps({cat: members for cat, members in ITP_CATEGORIES.items()})

# ══════════════════════════════════════════════════════════════════════════════
# ZIMS data
# ══════════════════════════════════════════════════════════════════════════════
zims = pd.read_csv(ZIMS_CSV)
try:
    tax = pd.read_csv(TAX_CSV)
    merge_cols = ["binSpecies","genus","family","order"] + (["common_name"] if "common_name" in tax.columns else [])
    zims = zims.merge(tax[merge_cols], on="binSpecies", how="left")
    zims["genus"]  = zims["genus"].fillna(zims["binSpecies"].str.split().str[0])
    zims["family"] = zims["family"].fillna("Unknown")
    zims["order"]  = zims["order"].fillna("Unknown")
    if "common_name" not in zims.columns: zims["common_name"] = ""
    zims["common_name"] = zims["common_name"].fillna("")
except:
    zims["genus"] = zims["binSpecies"].str.split().str[0]
    zims["family"] = zims["order"] = zims["common_name"] = "Unknown"

zims["lx"] = np.log10(zims["rho_beta"].clip(lower=1e-6))
zims["ly"] = np.log10(zims["rho_eta"].clip(lower=1e-6))
zims["lz"] = np.log10(zims["rho_eps"].clip(lower=1e-6))
s99 = zims["s"].clip(lower=1e-3).quantile(0.99)
zims["log_s"]   = np.log10(zims["s"].clip(lower=1e-3, upper=s99))
zims["log_L"]   = np.log10(zims["L"].clip(lower=1))
zims["log_bxe"] = np.log10(zims["beta_xc_eps"].clip(lower=1e-4))

CLASS_COLOR = {"Mammalia":"#e6194b","Aves":"#3cb44b","Reptilia":"#4363d8",
               "Amphibia":"#f58231","Chondrichthyes":"#911eb4"}
CLASS_EMOJI = {"Mammalia":"🦁","Aves":"🦅","Reptilia":"🦎",
               "Amphibia":"🐸","Chondrichthyes":"🦈"}
classes = sorted(zims["class"].unique())

zims_traces = []
for cls in classes:
    first = True
    for sex in ["f","m"]:
        for arm in ["removal","production"]:
            sub = zims[(zims["class"]==cls)&(zims["sex"]==sex)&(zims["arm"]==arm)]
            if sub.empty: continue
            hover = [
                f"<b>{r['binSpecies']}</b>"
                + (f" ({r['common_name']})" if r.get('common_name') else "")
                + f" ({'♀' if sex=='f' else '♂'})<br>"
                f"{r['class']} · {r.get('order','?')} · {r.get('family','?')}<br>"
                f"L={r['L']:.0f} d  rms={r['rms']:.4f}  arm={r['arm']}  ndims={int(r['ndims'])}<br>"
                f"s={r['s']:.3f}  β·Xc/ε={r['beta_xc_eps']:.3f}  κ={r['kappa']:.4f}<br>"
                f"n={int(r['n']):,}  n_dead={int(r['n_dead']):,}"
                for _,r in sub.iterrows()
            ]
            zims_traces.append({
                "type":"scatter3d",
                "name":f"{CLASS_EMOJI[cls]} {cls}",
                "x":sub["lx"].tolist(),"y":sub["ly"].tolist(),"z":sub["lz"].tolist(),
                "text":hover,"hoverinfo":"text","mode":"markers",
                "marker":{
                    "color":CLASS_COLOR[cls],
                    "symbol":"circle" if sex=="f" else "square",
                    "size":[4]*len(sub),"opacity":0.8,
                    "line":{"color":"rgba(0,0,0,0.2)","width":0.3},
                },
                "legendgroup":f"zims_{cls}","showlegend":first,
                "meta":{"class":cls,"sex":sex,"arm":arm,
                        "rms":sub["rms"].tolist(),"genus":sub["genus"].tolist(),
                        "order":sub["order"].tolist(),"family":sub["family"].tolist(),
                        "kappa":sub["kappa"].tolist(),
                        "class_color":CLASS_COLOR[cls],"dataset":"zims"},
            })
            first = False


# ── 2D data ───────────────────────────────────────────────────────────────────
traces_sL, traces_id = [], []
for cls in classes:
    sub = zims[zims["class"]==cls]
    hov = [f"<b>{r['binSpecies']}</b> ({'♀' if r['sex']=='f' else '♂'})<br>L={r['L']:.0f} d  s={r['s']:.3f}  rms={r['rms']:.4f}" for _,r in sub.iterrows()]
    traces_sL.append({"type":"scatter","name":f"{CLASS_EMOJI[cls]} {cls}",
                      "x":sub["log_L"].tolist(),"y":sub["log_s"].tolist(),
                      "mode":"markers","text":hov,"hoverinfo":"text",
                      "marker":{"color":CLASS_COLOR[cls],"size":5,"opacity":0.7},
                      "legendgroup":f"z2_{cls}","showlegend":True})
    traces_id.append({"type":"scatter","name":f"{CLASS_EMOJI[cls]} {cls}",
                      "x":sub["log_bxe"].tolist(),"y":sub["log_s"].tolist(),
                      "mode":"markers","hoverinfo":"skip",
                      "marker":{"color":CLASS_COLOR[cls],"size":5,"opacity":0.7},
                      "legendgroup":f"z2_{cls}","showlegend":False})
# ITP on sL
for sex in ["f","m"]:
    sub = itp[itp["sex"]==sex]
    traces_sL.append({"type":"scatter","name":f"ITP {'♀' if sex=='f' else '♂'}",
                      "x":np.log10(sub["L"]).tolist(),"y":np.log10(sub["s"].clip(lower=1e-3)).tolist(),
                      "mode":"markers","marker":{"color":SEX_COLOR[sex],"symbol":"star","size":8,"opacity":0.9,"line":{"color":"white","width":0.5}},
                      "text":[f"<b>{r['group']}</b> ITP ({'♀' if sex=='f' else '♂'})<br>L={r['L']:.0f} d  s={r['s']:.2f}" for _,r in sub.iterrows()],
                      "hoverinfo":"text","showlegend":True,"legendgroup":"itp2"})

# ── Phylo data ────────────────────────────────────────────────────────────────
coords  = zims[["lx","ly","lz"]].values.astype(np.float64)
n = len(zims)
i_idx, j_idx = np.triu_indices(n, k=1)
diff = coords[i_idx]-coords[j_idx]
dists_all = np.sqrt((diff**2).sum(axis=1))
genera  = zims["genus"].values
families= zims["family"].values
orders  = zims["order"].values
cls_arr = zims["class"].values
same_gen = genera[i_idx]==genera[j_idx]
same_fam = families[i_idx]==families[j_idx]
same_ord = orders[i_idx]==orders[j_idx]
same_cls = cls_arr[i_idx]==cls_arr[j_idx]
level = np.where(same_gen,0,np.where(same_fam,1,np.where(same_ord,2,np.where(same_cls,3,4))))
rng = np.random.default_rng(42)
phylo_traces = []
for lv,label,color in [(0,"Same genus","#ff3366"),(1,"Same family","#ff8800"),
                        (2,"Same order","#ffdd00"),(3,"Same class","#44cc44"),(4,"Diff class","#4499ff")]:
    mask = level==lv
    d = dists_all[mask]
    if not len(d): continue
    idx = rng.choice(len(d), min(8000,len(d)), replace=False)
    phylo_traces.append({"type":"violin","name":label,"y":d[idx].tolist(),
                          "box":{"visible":True},"meanline":{"visible":True},
                          "fillcolor":color,"opacity":0.7,"line":{"color":color},"points":False})

# ── Sunburst data ─────────────────────────────────────────────────────────────
def hex_rgba(h, a):
    r,g,b = int(h[1:3],16),int(h[3:5],16),int(h[5:7],16)
    return f"rgba({r},{g},{b},{a})"

sun_ids, sun_labels, sun_parents, sun_values, sun_colors = (
    ["root"],["All species"],[""],
    [int(zims["binSpecies"].nunique())],
    ["#374151"],
)
for cls in sorted(zims["class"].unique()):
    cc = CLASS_COLOR[cls]
    sub_c = zims[zims["class"]==cls]
    sun_ids.append(cls); sun_labels.append(cls); sun_parents.append("root")
    sun_values.append(int(sub_c["binSpecies"].nunique())); sun_colors.append(hex_rgba(cc,0.90))
    for ord_ in sorted(sub_c["order"].dropna().unique()):
        if not ord_ or ord_ in ("Unknown",""):  continue
        sub_o = sub_c[sub_c["order"]==ord_]
        oid = f"{cls}/{ord_}"
        sun_ids.append(oid); sun_labels.append(ord_); sun_parents.append(cls)
        sun_values.append(int(sub_o["binSpecies"].nunique())); sun_colors.append(hex_rgba(cc,0.72))
        for fam in sorted(sub_o["family"].dropna().unique()):
            if not fam or fam in ("Unknown",""): continue
            sub_f = sub_o[sub_o["family"]==fam]
            fid = f"{cls}/{ord_}/{fam}"
            sun_ids.append(fid); sun_labels.append(fam); sun_parents.append(oid)
            sun_values.append(int(sub_f["binSpecies"].nunique())); sun_colors.append(hex_rgba(cc,0.55))
            for gen in sorted(sub_f["genus"].dropna().unique()):
                sub_g = sub_f[sub_f["genus"]==gen]
                gid = f"{cls}/{ord_}/{fam}/{gen}"
                sun_ids.append(gid); sun_labels.append(gen); sun_parents.append(fid)
                sun_values.append(int(sub_g["binSpecies"].nunique())); sun_colors.append(hex_rgba(cc,0.40))

sunburst_trace = {
    "type":"sunburst",
    "ids":sun_ids,"labels":sun_labels,"parents":sun_parents,"values":sun_values,
    "branchvalues":"total","maxdepth":3,
    "marker":{"colors":sun_colors,"line":{"color":"#1f2937","width":0.5}},
    "textfont":{"color":"white","size":11},
    "insidetextfont":{"color":"white","size":11},
    "outsidetextfont":{"color":"white","size":11},
    "hovertemplate":"<b>%{label}</b><br>%{value} species<extra></extra>",
    "insidetextorientation":"radial",
    "leaf":{"opacity":0.85},
}
sunburst_json = json.dumps(sunburst_trace)

# ── Taxonomy hierarchy for cascading dropdowns ───────────────────────────────
tax_hier = {}
for cls in sorted(zims["class"].unique()):
    sub_c = zims[zims["class"]==cls]
    orders = {}
    for ord_ in sorted(sub_c.loc[~sub_c["order"].isin(["Unknown",""]), "order"].unique()):
        sub_o = sub_c[sub_c["order"]==ord_]
        families = {}
        for fam in sorted(sub_o.loc[~sub_o["family"].isin(["Unknown",""]), "family"].unique()):
            sub_f = sub_o[sub_o["family"]==fam]
            families[fam] = sorted(sub_f["genus"].unique().tolist())
        orders[ord_] = families
    tax_hier[cls] = orders
tax_hier_json = json.dumps(tax_hier)

# ── Collapsible taxonomy tree HTML ───────────────────────────────────────────
def _tcnt(zims, **kw):
    m = zims
    for col, val in kw.items():
        m = m[m[col]==val]
    return int(m["binSpecies"].nunique())

def tree_li(path, label, count, color, children_html, weight="normal"):
    tog = (f'<span class="tt" onclick="ttTog(this)">&#9656;</span>'
           if children_html else '<span class="tt-sp"></span>')
    cb  = f'<input type="checkbox" class="tcb" data-path="{path}" onchange="tcbChange(this)">'
    lbl = (f'<span class="tlbl" style="color:{color};font-weight:{weight};">{label}</span>')
    cnt = f'<span class="tcnt">{count}</span>'
    inner = (f'<ul class="tch">{children_html}</ul>' if children_html else "")
    return f'<li>{tog}{cb}{lbl}{cnt}{inner}</li>'

tree_html = ""
for cls in sorted(tax_hier.keys()):
    cc = CLASS_COLOR[cls]; em = CLASS_EMOJI[cls]
    orders_html = ""
    for ord_ in sorted(tax_hier[cls].keys()):
        fams_html = ""
        for fam in sorted(tax_hier[cls][ord_].keys()):
            genera = tax_hier[cls][ord_][fam]
            gens_html = ""
            for gen in sorted(genera):
                n = _tcnt(zims, **{"class":cls,"order":ord_,"family":fam,"genus":gen})
                gens_html += tree_li(
                    f"{cls}/{ord_}/{fam}/{gen}",
                    f"<i>{gen}</i>", n, "#d1d5db", "")
            n = _tcnt(zims, **{"class":cls,"order":ord_,"family":fam})
            fams_html += tree_li(
                f"{cls}/{ord_}/{fam}",
                f"{em}&thinsp;<i>{fam}</i>", n, "#e5e7eb", gens_html)
        n = _tcnt(zims, **{"class":cls,"order":ord_})
        orders_html += tree_li(
            f"{cls}/{ord_}",
            f"{em} {ord_}", n, "#f3f4f6", fams_html, weight="600")
    n = _tcnt(zims, **{"class":cls})
    tree_html += tree_li(
        cls, f"{em} {cls}", n, cc, orders_html, weight="700")

tree_html = f'<ul class="taxon-tree">{tree_html}</ul>'

# ── Naveh species checkboxes HTML ─────────────────────────────────────────────
naveh_cb_html = ""
for idx, name, em in NAVEH_SPECIES:
    naveh_cb_html += (f'<label class="cb-row">'
                      f'<input type="checkbox" class="naveh-cb" value="{idx}">'
                      f' {em} {name}</label>\n')

# ══════════════════════════════════════════════════════════════════════════════
# JSON payloads
# ══════════════════════════════════════════════════════════════════════════════
itp_traces_json   = json.dumps(itp_traces)
zims_traces_json  = json.dumps(zims_traces)
traces_sL_json    = json.dumps(traces_sL)
traces_id_json    = json.dumps(traces_id)
phylo_traces_json = json.dumps(phylo_traces)
naveh_idx_json    = json.dumps(NAVEH_IDX)
class_color_json  = json.dumps(CLASS_COLOR)

# ══════════════════════════════════════════════════════════════════════════════
# HTML injection
# ══════════════════════════════════════════════════════════════════════════════
inject = f"""
<style>
*{{box-sizing:border-box;}}
/* ── tab bar ── */
#sr-tabs{{
  position:fixed;top:0;left:0;right:0;height:42px;z-index:10000;
  background:#111827;border-bottom:2px solid #374151;
  display:flex;align-items:center;padding:0 10px;gap:4px;
}}
.sr-tab{{padding:8px 16px;cursor:pointer;background:transparent;border:none;
         color:#9ca3af;font-size:13px;font-weight:600;border-radius:6px;
         transition:all .15s;font-family:Arial,sans-serif;
         flex-shrink:0;white-space:nowrap;}}
.sr-tab:hover{{color:#e5e7eb;background:#374151;}}
.sr-tab.active{{color:#60a5fa;background:#1f2937;}}

/* ── tip bar (below tabs, always visible) ── */
#tip-bar{{
  position:fixed;top:42px;left:0;right:0;height:32px;z-index:9800;
  background:#0f1e3a;border-bottom:1px solid #1e3a5f;
  display:flex;align-items:center;padding:0 14px 0 284px;
  font-size:10.5px;color:#7eb8f7;font-family:Arial,sans-serif;
  white-space:nowrap;overflow:hidden;gap:0;
}}
#tip-bar .tsep{{color:#2d4a7a;margin:0 10px;}}

/* ── push main plot below tip bar; start invisible to hide load flash ── */
#{DIV_ID}{{margin-top:74px !important;height:calc(100vh - 74px) !important;
           opacity:0;transition:opacity .4s ease;}}

/* ── overlay views (Phylo) — start below tip bar ── */
.sr-view{{
  display:none;position:fixed;top:74px;left:270px;right:0;bottom:0;
  z-index:9000;background:#111827;
}}
.sr-view.active{{display:flex;}}
#sr-plot-sunburst,#sr-plot-phylo{{flex:1;min-height:0;}}
#sr-phylo-cap{{padding:7px 14px;font-size:10px;color:#6b7280;
               background:#1f2937;border-top:1px solid #374151;
               font-family:Arial,sans-serif;}}

/* ── filter panel (left) ── */
#sr-panel{{
  position:fixed;top:74px;left:0;bottom:22px;width:270px;z-index:9500;
  background:#1f2937;border-right:1px solid #374151;
  overflow-y:auto;padding:10px;font-size:12px;color:#e5e7eb;
  font-family:Arial,sans-serif;
}}

/* ── tutorial modal ── */
#tut-overlay{{
  display:none;position:fixed;inset:0;background:rgba(0,0,0,.75);
  z-index:99999;align-items:center;justify-content:center;
}}
#tut-box{{
  background:#1f2937;border:1px solid #374151;border-radius:12px;
  padding:26px 28px 20px;max-width:520px;width:92%;
  font-family:Arial,sans-serif;
}}
#tut-box h3{{margin:0 0 14px;color:#93c5fd;font-size:16px;}}
.tut-sec{{margin-bottom:14px;}}
.tut-sec b{{color:#fbbf24;font-size:11px;display:block;margin-bottom:3px;}}
.tut-sec p{{margin:0;font-size:11px;color:#d1d5db;line-height:1.7;}}
#tut-footer{{display:flex;align-items:center;justify-content:space-between;
             margin-top:18px;padding-top:14px;border-top:1px solid #374151;}}
.tut-dismiss{{font-size:11px;color:#6b7280;cursor:pointer;display:flex;align-items:center;gap:5px;}}
#tut-close{{padding:8px 22px;background:#2563eb;border:none;border-radius:6px;
            color:white;font-size:13px;font-weight:600;cursor:pointer;}}
#tut-close:hover{{background:#1d4ed8;}}
.help-btn{{
  margin-left:auto;padding:5px 10px;background:transparent;
  border:1px solid #374151;border-radius:5px;color:#6b7280;
  font-size:11px;cursor:pointer;font-family:Arial,sans-serif;
}}
.help-btn:hover{{color:#e5e7eb;border-color:#4b5563;background:#374151;}}

/* ── collapsible taxonomy tree ── */
ul.taxon-tree,ul.tch{{list-style:none;margin:0;padding:0;}}
ul.tch{{padding-left:14px;display:none;}}
ul.taxon-tree>li{{margin-bottom:1px;}}
ul.taxon-tree li{{line-height:1.6;}}
.tt{{display:inline-block;width:14px;cursor:pointer;font-size:9px;
     color:#6b7280;user-select:none;text-align:center;transition:transform .15s;}}
.tt-sp{{display:inline-block;width:14px;}}
.tcb{{margin:0 3px 0 0;cursor:pointer;flex-shrink:0;vertical-align:middle;}}
.tlbl{{cursor:default;font-size:11px;vertical-align:middle;}}
.tcnt{{font-size:9px;color:#6b7280;margin-left:3px;vertical-align:middle;}}
.sec{{margin-bottom:2px;border-radius:6px;overflow:hidden;background:#111827;}}
.sec-hdr{{
  display:flex;justify-content:space-between;align-items:center;
  padding:8px 10px;cursor:pointer;font-weight:700;font-size:12px;
  color:#93c5fd;background:#1e3a5f;border-radius:6px;
  user-select:none;transition:background .15s;
}}
.sec-hdr:hover{{background:#1a4a7a;}}
.sec-hdr .arrow{{font-size:10px;color:#6b7280;transition:transform .2s;}}
.sec-hdr.open .arrow{{transform:rotate(90deg);}}
.sec-body{{padding:8px;display:none;}}
.sec-body.open{{display:block;}}
h5{{margin:6px 0 3px;font-size:10px;color:#6b7280;text-transform:uppercase;letter-spacing:.5px;font-weight:700;}}
.cb-row{{display:flex;align-items:center;gap:5px;padding:2px 4px;margin-bottom:1px;
         border-radius:4px;cursor:pointer;font-size:11px;}}
.cb-row:hover{{background:#374151;}}
.cb-row input{{margin:0;cursor:pointer;flex-shrink:0;}}
.row2{{display:grid;grid-template-columns:1fr 1fr;gap:3px;}}
.cb-grid{{display:grid;grid-template-columns:1fr 1fr;gap:2px;}}
select{{width:100%;padding:4px;background:#374151;border:1px solid #4b5563;
        border-radius:5px;color:#e5e7eb;font-size:11px;margin-bottom:4px;}}
input[type=range]{{width:100%;accent-color:#60a5fa;margin:2px 0;}}
.rms-val{{color:#fbbf24;font-size:10px;float:right;}}
.btn-row{{display:flex;gap:5px;margin-bottom:6px;}}
.mini-btn{{flex:1;padding:3px 0;background:#374151;border:1px solid #4b5563;
           border-radius:4px;color:#e5e7eb;font-size:10px;cursor:pointer;}}
.mini-btn:hover{{background:#4b5563;}}
.note{{font-size:10px;color:#6b7280;line-height:1.6;margin-top:6px;}}
.good{{color:#34d399;}} .warn{{color:#fbbf24;}} .bad{{color:#f87171;}}
</style>

<!-- Tab bar -->
<div id="sr-tabs">
  <button class="sr-tab active" onclick="srTab('3d')">3D Manifold</button>
  <button class="sr-tab"        onclick="srTab('phylo')">Phylo Tree</button>
  <button class="help-btn"      onclick="openTutorial()">? Help</button>
</div>

<!-- Tip bar -->
<div id="tip-bar">
  <span>🖱 Drag: rotate &nbsp;·&nbsp; scroll: zoom</span>
  <span class="tsep">|</span>
  <span>Left panel: check taxa, toggle data layers</span>
  <span class="tsep">|</span>
  <span>🌳 <b style="color:#93c5fd;">Phylo Tree</b> tab: click a sunburst wedge to focus the manifold on any taxon</span>
  <span class="tsep">|</span>
  <span>Multi-taxon: check multiple branches in the panel tree</span>
</div>

<!-- Tutorial modal -->
<div id="tut-overlay" onclick="if(event.target===this)closeTutorial()">
  <div id="tut-box">
    <h3>SR Manifold Viewer — Quick Start</h3>
    <div class="tut-sec">
      <b>Navigate the 3D scene</b>
      <p>Drag to rotate &nbsp;·&nbsp; scroll to zoom &nbsp;·&nbsp; double-click to reset camera.<br>
         Use <em>Show all / Hide all data</em> (top-right of plot) to toggle every data layer at once.</p>
    </div>
    <div class="tut-sec">
      <b>Show / hide data layers</b>
      <p>Left panel → <em>Manifold</em>: toggle the surface, ridge, and ω labels independently.<br>
         <em>Other species</em> (Naveh fits) and <em>NIA ITP</em> (mouse interventions): use their section toggles.<br>
         <em>ZIMS Zoo Animals</em>: check items in the taxon tree to show those species on the manifold.</p>
    </div>
    <div class="tut-sec">
      <b>Filter zoo animals by taxonomy</b>
      <p>Expand the tree in the ZIMS panel (▶ triangles open branches).<br>
         Check a class, order, family, or genus to display only those species.<br>
         Multiple checks = union — mix any taxa for a custom view.<br>
         <em>Check all</em> shows every species; <em>Uncheck all</em> hides them all.</p>
    </div>
    <div class="tut-sec">
      <b>Phylo Tree tab — one-click taxon focus</b>
      <p>Switch to the <em>Phylo Tree</em> tab and click any wedge of the sunburst chart.<br>
         That taxon is instantly checked in the left panel and the 3D manifold updates.<br>
         Click the centre ring (or <em>✕ clear filter</em>) to go back to the full view.<br>
         The violin chart (right) shows phylo-distance vs. proximity in SR parameter space.</p>
    </div>
    <div class="tut-sec">
      <b>Color mode</b>
      <p>Left panel → <em>Color</em>: switch between class/group coloring and κ (saturation) coloring for ITP &amp; ZIMS points.</p>
    </div>
    <div id="tut-footer">
      <label class="tut-dismiss">
        <input type="checkbox" id="tut-no-show"> Don't show again
      </label>
      <button id="tut-close" onclick="closeTutorial()">Got it!</button>
    </div>
  </div>
</div>

<!-- Phylo overlay -->
<div id="sr-view-phylo" class="sr-view" style="flex-direction:row;">
  <!-- sunburst panel -->
  <div style="flex:3;display:flex;flex-direction:column;min-width:0;border-right:1px solid #374151;">
    <div style="padding:6px 10px;background:#1f2937;flex-shrink:0;font-family:Arial,sans-serif;">
      <div style="display:flex;align-items:center;gap:8px;">
        <span style="color:#93c5fd;font-size:11px;font-weight:600;">Phylogenetic sunburst</span>
        <button id="taxon-reset" onclick="filterByTaxon('root')"
                style="display:none;padding:2px 8px;background:#374151;border:1px solid #4b5563;
                       border-radius:4px;color:#fbbf24;font-size:10px;cursor:pointer;flex-shrink:0;">
          ✕ clear filter</button>
      </div>
      <div style="font-size:10px;color:#6b7280;line-height:1.7;margin-top:3px;">
        <b style="color:#9ca3af;">How to use:</b>
        Click a wedge to check that taxon in the left panel and filter the 3D manifold to those species.
        Click the centre ring to clear. Scroll to zoom in/out. The filter also works independently
        from the left panel tree — check multiple branches there for a custom multi-taxon view.
      </div>
    </div>
    <div id="sr-plot-sunburst" style="flex:1;min-height:0;"></div>
  </div>
  <!-- violin panel -->
  <div style="flex:2;display:flex;flex-direction:column;min-width:0;">
    <div id="sr-plot-phylo" style="flex:1;min-height:0;"></div>
    <div id="sr-phylo-cap">
      Distance in log(ρ<sub>β</sub>,ρ<sub>η</sub>,ρ<sub>ε</sub>) vs taxonomic proximity.
      Violin = distribution · box = IQR · line = mean · ≤8 000 pairs/level.
    </div>
  </div>
</div>

<!-- Unified filter panel -->
<div id="sr-panel">

  <!-- 0. Color mode -->
  <div class="sec">
    <div class="sec-hdr" onclick="toggleSec(this)">
      🎨 Color <span class="arrow">▶</span>
    </div>
    <div class="sec-body">
      <h5>Color points by</h5>
      <select id="color-by" onchange="applyColorMode()">
        <option value="class">Class / group</option>
        <option value="kappa">κ (saturation)</option>
      </select>
      <div class="note" id="kappa-note" style="display:none;">
        log<sub>10</sub> κ · applies to ITP &amp; ZIMS points
      </div>
    </div>
  </div>

  <!-- 1. Manifold elements -->
  <div class="sec">
    <div class="sec-hdr" onclick="toggleSec(this)">
      🗺 Manifold <span class="arrow">▶</span>
    </div>
    <div class="sec-body">
      <label class="cb-row"><input type="checkbox" id="man-surface" checked onchange="applyManifold()"> Surface</label>
      <label class="cb-row"><input type="checkbox" id="man-ridge"   checked onchange="applyManifold()"> Ridge</label>
      <label class="cb-row"><input type="checkbox" id="man-omega"   checked onchange="applyManifold()"> ω labels</label>
      <div class="note" style="margin-top:6px;padding-top:6px;border-top:1px solid #1e3a5f;">
        <span style="color:#fbbf24;">⚠</span>
        Surface: κ=0 &nbsp;·&nbsp; Naveh fits: κ=0.5<br>
        ITP &amp; ZIMS: κ free, X<sub>c</sub>=1
      </div>
    </div>
  </div>

  <!-- 2. Naveh species -->
  <div class="sec">
    <div class="sec-hdr" onclick="toggleSec(this)">
      🔬 Other species <span class="arrow">▶</span>
    </div>
    <div class="sec-body">
      <div class="btn-row">
        <button class="mini-btn" onclick="setNaveh(true)">Show all</button>
        <button class="mini-btn" onclick="setNaveh(false)">Hide all</button>
      </div>
      <div class="cb-grid">
        {naveh_cb_html}
      </div>
    </div>
  </div>

  <!-- 3. ITP -->
  <div class="sec">
    <div class="sec-hdr" onclick="toggleSec(this)">
      🐭 NIA ITP <span class="arrow">▶</span>
    </div>
    <div class="sec-body">
      <button class="mini-btn" id="itp-toggle" onclick="toggleAllItp()"
              style="width:100%;margin-bottom:6px;background:#1e3a5f;border-color:#2563eb;color:#93c5fd;">
        Show all ITP</button>
      <h5>Intervention</h5>
      <select id="itp-grp" onchange="applyFilter()">{itp_opts}</select>
      <h5>Sex</h5>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="itp-f" checked onchange="applyFilter()"> <span style="color:#e8634a;">♀ Female</span></label>
        <label class="cb-row"><input type="checkbox" id="itp-m" checked onchange="applyFilter()"> <span style="color:#2c7bb6;">♂ Male</span></label>
      </div>
      <div class="note">● removal &nbsp; ◆ production<br>Larger/outlined = Control</div>
    </div>
  </div>

  <!-- 4. ZIMS -->
  <div class="sec">
    <div class="sec-hdr" onclick="toggleSec(this)">
      🌍 ZIMS Zoo Animals <span class="arrow">▶</span>
    </div>
    <div class="sec-body">
      <h5>Taxon (check to show)</h5>
      <div class="btn-row">
        <button class="mini-btn" onclick="setAllTree(true)">Check all</button>
        <button class="mini-btn" onclick="setAllTree(false)">Uncheck all</button>
      </div>
      <div style="max-height:260px;overflow-y:auto;border:1px solid #374151;border-radius:4px;padding:4px;background:#111827;">
        {tree_html}
      </div>
      <h5 style="margin-top:5px;">Sex</h5>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="zims-f" checked onchange="applyFilter()"> ♀ Female</label>
        <label class="cb-row"><input type="checkbox" id="zims-m" checked onchange="applyFilter()"> ♂ Male</label>
      </div>
      <h5 style="margin-top:5px;">Arm</h5>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="zims-rem" checked onchange="applyFilter()"> removal</label>
        <label class="cb-row"><input type="checkbox" id="zims-pro" checked onchange="applyFilter()"> production</label>
      </div>
      <h5 style="margin-top:5px;">Max RMS <span class="rms-val" id="rms-disp">all</span></h5>
      <input type="range" id="rms-slider" min="0" max="160" value="160" step="1"
             oninput="document.getElementById('rms-disp').textContent=this.value>=160?'all':(this.value/1000).toFixed(3);applyFilter()">
      <h5 style="margin-top:5px;">Dot size <span class="rms-val" id="sz-disp">4</span></h5>
      <input type="range" id="sz-slider" min="1" max="12" value="4" step="1"
             oninput="document.getElementById('sz-disp').textContent=this.value;applyFilter()">
      <div class="note">
        ● ♀ circle &nbsp; ■ ♂ square<br>
        Median RMS: <span class="good">{zims['rms'].median():.4f}</span> &nbsp;
        ndims=5: <span class="warn">{int((zims['ndims']==5).sum())}</span> &nbsp;
        Poor fits: <span class="bad">{int((zims['rms']>0.10).sum())}</span>
      </div>
    </div>
  </div>

</div><!-- /sr-panel -->

<script>
(function(){{
var DIVID         = "{DIV_ID}";
var gd            = document.getElementById(DIVID);
var itpTraces     = {itp_traces_json};
var zimsTraces    = {zims_traces_json};
var phyloTr       = {phylo_traces_json};
var catMap        = {itp_cat_map_json};
var navehIdx      = {naveh_idx_json};
var clsColor      = {class_color_json};
var sunburstTrace = {sunburst_json};
var taxHier       = {tax_hier_json};
var itpStart, zimsStart;
var tabPhyloDone=false;

// ── tutorial modal ────────────────────────────────────────────────────────────
window.openTutorial = function() {{
  document.getElementById('tut-overlay').style.display = 'flex';
}};
window.closeTutorial = function() {{
  document.getElementById('tut-overlay').style.display = 'none';
  if (document.getElementById('tut-no-show').checked) {{
    try {{ localStorage.setItem('sr_tut_seen','1'); }} catch(e) {{}}
  }}
}};
(function() {{
  var seen = false;
  try {{ seen = localStorage.getItem('sr_tut_seen')==='1'; }} catch(e) {{}}
  if (!seen) openTutorial();
}})();

// ── section collapse ──────────────────────────────────────────────────────────
window.toggleSec = function(hdr) {{
  hdr.classList.toggle('open');
  hdr.nextElementSibling.classList.toggle('open');
}};

// ── taxonomy tree ──────────────────────────────────────────────────────────────
window.ttTog = function(span) {{
  var li = span.parentElement;
  var ch = li.querySelector('ul.tch');
  if (!ch) return;
  var open = ch.style.display === 'block';
  ch.style.display = open ? 'none' : 'block';
  span.innerHTML = open ? '&#9656;' : '&#9662;';
}};
window.tcbChange = function(cb) {{
  // cascade to all descendants within same <li>
  cb.parentElement.querySelectorAll('.tcb').forEach(function(d){{ d.checked=cb.checked; }});
  applyFilter();
}};
window.setAllTree = function(v) {{
  document.querySelectorAll('.tcb').forEach(function(cb){{ cb.checked=v; }});
  applyFilter();
}};

// ── tab switching ─────────────────────────────────────────────────────────────
window.srTab = function(tab) {{
  document.querySelectorAll('.sr-tab').forEach(function(b,i){{
    b.classList.toggle('active', ['3d','phylo'][i]===tab);
  }});
  document.getElementById('sr-view-phylo').classList.toggle('active', tab==='phylo');
  if (tab==='phylo' && !tabPhyloDone){{ tabPhyloDone=true;
    requestAnimationFrame(function(){{requestAnimationFrame(initPhylo);}}); }}
  if (tab==='3d') Plotly.Plots.resize(DIVID);
}};

function initPhylo() {{
  // Set explicit heights so Plotly can measure the containers
  var avail = window.innerHeight - 74;
  var sbEl = document.getElementById('sr-plot-sunburst');
  var viEl = document.getElementById('sr-plot-phylo');
  sbEl.style.height = Math.round(avail * 0.97) + 'px';
  viEl.style.height = Math.round(avail * 0.87) + 'px';
  // Sunburst
  Plotly.newPlot('sr-plot-sunburst', [sunburstTrace], {{
    paper_bgcolor:'#111827',
    font:{{color:'#e5e7eb',family:'Arial,sans-serif',size:11}},
    margin:{{l:0,r:0,t:0,b:0}},
  }}, {{responsive:true}});
  document.getElementById('sr-plot-sunburst').on('plotly_click', function(d){{
    if (!d.points || !d.points.length) return;
    var id = d.points[0].id;
    filterByTaxon(id);
  }});
  // Violin
  Plotly.newPlot('sr-plot-phylo', phyloTr, {{
    paper_bgcolor:'#111827',plot_bgcolor:'#1f2937',
    font:{{color:'#e5e7eb',family:'Arial,sans-serif'}},
    yaxis:{{title:'Distance in log ρ-space',gridcolor:'#374151',color:'#9ca3af'}},
    xaxis:{{color:'#9ca3af',tickfont:{{size:9}}}},
    title:{{text:'Distance vs Taxonomic Proximity',font:{{color:'#93c5fd',size:12}}}},
    violinmode:'overlay',showlegend:false,
    margin:{{l:55,r:10,t:40,b:60}},
  }},{{responsive:true}});
}}

// ── taxonomy filter (from sunburst click → drives tree checkboxes) ───────────
window.filterByTaxon = function(id) {{
  var btn = document.getElementById('taxon-reset');
  // Uncheck everything first
  document.querySelectorAll('.tcb').forEach(function(cb){{ cb.checked=false; }});
  if (!id || id==='root') {{
    btn.style.display = 'none';
  }} else {{
    // Check just this node + expand its ancestors
    var cb = document.querySelector('.tcb[data-path="'+id+'"]');
    if (cb) {{
      cb.checked = true;
      // Also cascade to descendants
      cb.parentElement.querySelectorAll('.tcb').forEach(function(d){{ d.checked=true; }});
      // Expand the panel section so the tree is visible
      var sec = document.querySelector('#sr-panel .sec:last-child .sec-hdr');
      if (sec && !sec.classList.contains('open')) toggleSec(sec);
    }}
    btn.style.display = 'inline-block';
  }}
  applyFilter();
}};

// ── inject ITP + ZIMS traces ──────────────────────────────────────────────────
function tryAdd() {{
  if (!gd._fullLayout) {{ setTimeout(tryAdd, 200); return; }}
  itpStart = gd.data.length;
  Plotly.addTraces(DIVID, itpTraces).then(function() {{
    zimsStart = gd.data.length;
    Plotly.addTraces(DIVID, zimsTraces).then(function() {{
      // extend show/hide buttons to cover all new traces
      var allNew=[];
      for (var i=itpStart; i<gd.data.length; i++) allNew.push(i);
      var ext = navehIdx.concat(allNew);
      Plotly.relayout(DIVID, {{
        showlegend: false,
        'updatemenus[0].buttons[0].args': [{{visible:true}},        ext],
        'updatemenus[0].buttons[1].args': [{{visible:'legendonly'}}, ext],
      }});
      applyNaveh();
      applyFilter();
      applyColorMode();
      document.getElementById(DIVID).style.opacity = '1';
      // Hook into the original "show all / hide all data" Plotly updatemenu buttons
      gd.on('plotly_buttonclicked', function(e) {{
        if (e.button && e.button.args && e.button.args[0] &&
            e.button.args[0].visible === true) {{
          // "Show all" button — reset our filters so everything is visible
          requestAnimationFrame(function() {{
            setAllTree(true);     // check all tree nodes → all ZIMS visible
            setNaveh(true);       // show all Naveh species
            if (!itpVisible) {{
              itpVisible = true;
              document.getElementById('itp-toggle').textContent = 'Hide all ITP';
            }}
            applyFilter();
          }});
        }}
      }});
    }});
  }});
}}
tryAdd();

// ── manifold elements filter ──────────────────────────────────────────────────
window.applyManifold = function() {{
  Plotly.restyle(DIVID, {{visible: [
    document.getElementById('man-surface').checked,
    document.getElementById('man-ridge').checked,
    document.getElementById('man-omega').checked,
  ]}}, [0, 1, 2]);
}};

// ── Naveh species filter ──────────────────────────────────────────────────────
window.setNaveh = function(show) {{
  document.querySelectorAll('.naveh-cb').forEach(function(cb){{ cb.checked=show; }});
  applyNaveh();
}};
function applyNaveh() {{
  var vis=[], idx=[];
  document.querySelectorAll('.naveh-cb').forEach(function(cb){{
    idx.push(parseInt(cb.value));
    vis.push(cb.checked);
  }});
  Plotly.restyle(DIVID, {{visible: vis}}, idx);
}}
document.querySelectorAll('.naveh-cb').forEach(function(cb){{
  cb.addEventListener('change', applyNaveh);
}});

// ── κ color mode ─────────────────────────────────────────────────────────────
var _kappaMin, _kappaMax;
function kappaRange() {{
  if (_kappaMin !== undefined) return;
  var vals = [];
  zimsTraces.forEach(function(t) {{
    t.meta.kappa.forEach(function(k) {{ if (k>0) vals.push(Math.log10(k)); }});
  }});
  itpTraces.forEach(function(t) {{
    (t.meta.kappa||[]).forEach(function(k) {{ if (k>0) vals.push(Math.log10(k)); }});
  }});
  _kappaMin = Math.min.apply(null,vals);
  _kappaMax = Math.max.apply(null,vals);
}}
window.applyColorMode = function() {{
  if (zimsStart===undefined) return;
  var mode = document.getElementById('color-by').value;
  document.getElementById('kappa-note').style.display = mode==='kappa' ? 'block':'none';
  var nZ = zimsTraces.length, nI = itpTraces.length;
  var zIdx = zimsTraces.map(function(_,i){{ return zimsStart+i; }});
  var iIdx = itpTraces.map(function(_,i){{ return itpStart+i; }});
  if (mode==='kappa') {{
    kappaRange();
    var zColors = zimsTraces.map(function(t){{ return t.meta.kappa.map(function(k){{ return k>0?Math.log10(k):_kappaMin; }}); }});
    var iColors = itpTraces.map(function(t){{ return (t.meta.kappa||[]).map(function(k){{ return k>0?Math.log10(k):_kappaMin; }}); }});
    var zShow = zimsTraces.map(function(_,i){{ return i===0; }});
    var iShow = itpTraces.map(function(){{ return false; }});
    Plotly.restyle(DIVID, {{
      'marker.color': zColors,
      'marker.colorscale': Array(nZ).fill('Plasma'),
      'marker.cmin': Array(nZ).fill(_kappaMin),
      'marker.cmax': Array(nZ).fill(_kappaMax),
      'marker.showscale': zShow,
      'marker.colorbar': [{{title:{{text:'log₁₀ κ'}},thickness:14,len:0.6,x:0.82,y:0.5,yanchor:'middle',bgcolor:'rgba(31,41,55,0.85)',outlinecolor:'#4b5563',tickfont:{{color:'#e5e7eb',size:10}},titlefont:{{color:'#93c5fd',size:11}}}}].concat(Array(nZ-1).fill(null)),
    }}, zIdx);
    Plotly.restyle(DIVID, {{
      'marker.color': iColors,
      'marker.colorscale': Array(nI).fill('Plasma'),
      'marker.cmin': Array(nI).fill(_kappaMin),
      'marker.cmax': Array(nI).fill(_kappaMax),
      'marker.showscale': iShow,
    }}, iIdx);
  }} else {{
    var zColors = zimsTraces.map(function(t){{ return t.meta.class_color; }});
    var iColors = itpTraces.map(function(t){{ return t.meta.sex_color; }});
    Plotly.restyle(DIVID, {{
      'marker.color': zColors,
      'marker.colorscale': Array(nZ).fill(null),
      'marker.showscale': Array(nZ).fill(false),
    }}, zIdx);
    Plotly.restyle(DIVID, {{
      'marker.color': iColors,
      'marker.colorscale': Array(nI).fill(null),
      'marker.showscale': Array(nI).fill(false),
    }}, iIdx);
  }}
}};

// ── ITP show/hide toggle ─────────────────────────────────────────────────────
var itpVisible = false;
window.toggleAllItp = function() {{
  itpVisible = !itpVisible;
  document.getElementById('itp-toggle').textContent = itpVisible ? 'Hide all ITP' : 'Show all ITP';
  applyFilter();
}};

// ── ITP + ZIMS filter ─────────────────────────────────────────────────────────
function itpGroupsFor(val) {{
  if (val==='all') return null;
  if (val.startsWith('cat:')) return catMap[val.slice(4)]||[];
  return [val];
}}

window.applyFilter = function() {{
  if (itpStart===undefined || zimsStart===undefined) return;

  // ITP
  var itpAllowed = itpGroupsFor(document.getElementById('itp-grp').value);
  var itpF = document.getElementById('itp-f').checked;
  var itpM = document.getElementById('itp-m').checked;
  var itpVis = itpTraces.map(function(t) {{
    if (!itpVisible) return false;
    var m=t.meta;
    return ((m.sex==='f'&&itpF)||(m.sex==='m'&&itpM)) &&
           (itpAllowed===null || itpAllowed.indexOf(m.group)!==-1);
  }});
  Plotly.restyle(DIVID, {{visible:itpVis}},
    itpTraces.map(function(_,i){{ return itpStart+i; }}));

  // ZIMS — tree-based filter.  Empty checkedSet = show nothing (empty default).
  var checkedSet = new Set();
  document.querySelectorAll('.tcb:checked').forEach(function(cb){{
    checkedSet.add(cb.dataset.path);
  }});

  var zF  = document.getElementById('zims-f').checked;
  var zM  = document.getElementById('zims-m').checked;
  var zR  = document.getElementById('zims-rem').checked;
  var zP  = document.getElementById('zims-pro').checked;
  var rms = document.getElementById('rms-slider').value/1000;
  var rmsAll = rms >= 0.159;
  var dotSz = parseInt(document.getElementById('sz-slider').value)||4;

  var zVis=[], zSz=[], zIdx=[];
  zimsTraces.forEach(function(t,i){{
    var m=t.meta;
    var clsPrefix = m.class+'/';
    var traceVis = false;
    checkedSet.forEach(function(p){{ if (p===m.class || p.indexOf(clsPrefix)===0) traceVis=true; }});
    zVis.push(traceVis &&
              ((m.sex==='f'&&zF)||(m.sex==='m'&&zM)) &&
              ((m.arm==='removal'&&zR)||(m.arm==='production'&&zP)));
    zSz.push(m.rms.map(function(r,j){{
      var cls=m.class, ord=m.order[j], fam=m.family[j], gen=m.genus[j];
      var pass = (rmsAll||r<=rms) && (
        checkedSet.has(cls) ||
        checkedSet.has(cls+'/'+ord) ||
        checkedSet.has(cls+'/'+ord+'/'+fam) ||
        checkedSet.has(cls+'/'+ord+'/'+fam+'/'+gen)
      );
      return pass ? dotSz : 0;
    }}));
    zIdx.push(zimsStart+i);
  }});
  Plotly.restyle(DIVID, {{visible:zVis}}, zIdx);
  Plotly.restyle(DIVID, {{'marker.size':zSz}}, zIdx);
}};

document.getElementById('itp-grp').addEventListener('change', applyFilter);
}})();
</script>
"""

with open(HTML_IN) as f:
    html = f.read()
html = html.replace("</body>", inject + "\n</body>", 1)
with open(HTML_OUT, "w") as f:
    f.write(html)

import os
print(f"Done")
print(f"  ITP traces : {len(itp_traces)}")
print(f"  ZIMS traces: {len(zims_traces)}")
print(f"  File size  : {os.path.getsize(HTML_OUT)/1e6:.1f} MB")
print(f"  Saved: {HTML_OUT}")
