"""
Build the combined SR manifold viewer.
Reads index.html (clean manifold) and injects:
  - Tab bar: 3D Manifold | 2D Analysis | Phylo Analysis
  - Unified filter panel (ITP + ZIMS)
  - ITP traces
  - ZIMS traces
  - 2D s-vs-L and identifiable-space Plotly divs
  - Phylo violin Plotly div
Output: saturated-removal-manifold/viewer.html
"""
import pandas as pd, numpy as np, json, math

ITP_CSV  = (
    "/Users/Avimayo/Library/CloudStorage/GoogleDrive-mayoavi@gmail.com"
    "/.shortcut-targets-by-id/1JUcTeGwLojxa9QVWv-SdPTq9kAz4r-xR"
    "/SysMedBook2021/Chapter 6_Aging/ITP survival analysis/itp_fp_fits.csv"
)
ZIMS_CSV = "/Users/Avimayo/Documents/saturated_removal_paper/zims_all.csv"
TAX_CSV  = "/tmp/zims_taxonomy.csv"
HTML_IN  = "/Users/Avimayo/Documents/saturated-removal-manifold/index.html"
HTML_OUT = "/Users/Avimayo/Documents/saturated-removal-manifold/viewer.html"
DIV_ID   = "a1d89836-0daf-4582-ae54-c7c4aca77de6"

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
        symbols = [ARM_SYMBOL[a] for a in sub["arm"]]
        hover = [
            f"<b>{r['group']}</b> | {r['cohort']} | {'♀' if r['sex']=='f' else '♂'}<br>"
            f"L={r['L']:.0f} d  s={r['s']:.2f}  rms={r['rms']:.3f}<br>"
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
                "color":color,"symbol":symbols,
                "size":7 if is_ctrl else 5,"opacity":0.95,
                "line":{"color":"black" if is_ctrl else color,
                        "width":1.5 if is_ctrl else 0},
            },
            "legendgroup":grp,"showlegend":True,
            "meta":{"group":grp,"sex":sex,"dataset":"itp"},
        })

# ITP dropdown HTML
itp_opts_html = '<option value="all">— All interventions —</option>\n'
for cat, members in ITP_CATEGORIES.items():
    itp_opts_html += f'<optgroup label="── {cat} ──">\n'
    itp_opts_html += f'  <option value="cat:{cat}">★ All {cat}</option>\n'
    for g in members:
        itp_opts_html += f'  <option value="{g}">&nbsp;&nbsp;{g}</option>\n'
    itp_opts_html += '</optgroup>\n'

itp_cat_map_js = json.dumps({cat: members for cat, members in ITP_CATEGORIES.items()})

# ══════════════════════════════════════════════════════════════════════════════
# ZIMS data
# ══════════════════════════════════════════════════════════════════════════════
zims = pd.read_csv(ZIMS_CSV)
try:
    tax = pd.read_csv(TAX_CSV)
    zims = zims.merge(tax[["binSpecies","genus","family","order"]], on="binSpecies", how="left")
    zims["genus"]  = zims["genus"].fillna(zims["binSpecies"].str.split().str[0])
    zims["family"] = zims["family"].fillna("Unknown")
    zims["order"]  = zims["order"].fillna("Unknown")
except:
    zims["genus"]  = zims["binSpecies"].str.split().str[0]
    zims["family"] = "Unknown"; zims["order"] = "Unknown"

zims["lx"] = np.log10(zims["rho_beta"].clip(lower=1e-6))
zims["ly"] = np.log10(zims["rho_eta"].clip(lower=1e-6))
zims["lz"] = np.log10(zims["rho_eps"].clip(lower=1e-6))
s99 = zims["s"].clip(lower=1e-3).quantile(0.99)
zims["log_s"]   = np.log10(zims["s"].clip(lower=1e-3, upper=s99))
zims["log_L"]   = np.log10(zims["L"].clip(lower=1))
zims["log_bxe"] = np.log10(zims["beta_xc_eps"].clip(lower=1e-4))

CLASS_COLOR = {
    "Mammalia":      "#e6194b",
    "Aves":          "#3cb44b",
    "Reptilia":      "#4363d8",
    "Amphibia":      "#f58231",
    "Chondrichthyes":"#911eb4",
}
CLASS_EMOJI = {
    "Mammalia":"🦁","Aves":"🦅","Reptilia":"🦎","Amphibia":"🐸","Chondrichthyes":"🦈",
}
classes = sorted(zims["class"].unique())

# 3D ZIMS traces: one per (class, sex, arm) — class color never mutates
zims_traces = []
for cls in classes:
    first = True
    for sex in ["f","m"]:
        sx = "Female" if sex=="f" else "Male"
        for arm in ["removal","production"]:
            sub = zims[(zims["class"]==cls)&(zims["sex"]==sex)&(zims["arm"]==arm)]
            if sub.empty: continue
            hover = [
                f"<b>{r['binSpecies']}</b> ({'♀' if sex=='f' else '♂'})<br>"
                f"{r['class']} · {r.get('order','?')} · {r.get('family','?')}<br>"
                f"L={r['L']:.0f} d  rms={r['rms']:.4f}  arm={r['arm']}  ndims={int(r['ndims'])}<br>"
                f"s={r['s']:.3f}  β·Xc/ε={r['beta_xc_eps']:.3f}<br>"
                f"n={int(r['n']):,}  n_dead={int(r['n_dead']):,}"
                for _,r in sub.iterrows()
            ]
            symbol = "circle" if sex=="f" else "square"
            zims_traces.append({
                "type":"scatter3d",
                "name":f"{CLASS_EMOJI[cls]} {cls}",
                "x":sub["lx"].tolist(),"y":sub["ly"].tolist(),"z":sub["lz"].tolist(),
                "text":hover,"hoverinfo":"text","mode":"markers",
                "marker":{
                    "color":CLASS_COLOR[cls],"symbol":symbol,
                    "size":[3]*len(sub),"opacity":0.8,
                    "line":{"color":"rgba(0,0,0,0.2)","width":0.3},
                },
                "legendgroup":f"zims_{cls}","showlegend":first,
                "meta":{"class":cls,"sex":sex,"arm":arm,
                        "rms":sub["rms"].tolist(),
                        "genus":sub["genus"].tolist(),
                        "dataset":"zims"},
            })
            first = False

# Emoji centroid labels
for cls in classes:
    sub = zims[zims["class"]==cls]
    zims_traces.append({
        "type":"scatter3d","name":f"_{cls}_lbl",
        "x":[float(sub["lx"].mean())],"y":[float(sub["ly"].mean())],"z":[float(sub["lz"].mean())],
        "mode":"text","text":[CLASS_EMOJI[cls]],
        "textfont":{"size":20,"color":CLASS_COLOR[cls]},
        "hoverinfo":"skip","showlegend":False,
        "meta":{"class":cls,"sex":"lbl","arm":"lbl","rms":[0],"genus":[""],"dataset":"zims"},
    })

# ── 2D traces ──────────────────────────────────────────────────────────────
traces_sL, traces_id = [], []
for cls in classes:
    sub = zims[zims["class"]==cls]
    hov_sL = [f"<b>{r['binSpecies']}</b> ({'♀' if r['sex']=='f' else '♂'})<br>L={r['L']:.0f} d  s={r['s']:.3f}  rms={r['rms']:.4f}" for _,r in sub.iterrows()]
    hov_id = [f"<b>{r['binSpecies']}</b><br>s={r['s']:.3f}  β·Xc/ε={r['beta_xc_eps']:.3f}" for _,r in sub.iterrows()]
    base = {"mode":"markers","marker":{"color":CLASS_COLOR[cls],"size":5,"opacity":0.7},"legendgroup":f"zims2_{cls}"}
    traces_sL.append({**base,"type":"scatter","name":f"{CLASS_EMOJI[cls]} {cls}",
                      "x":sub["log_L"].tolist(),"y":sub["log_s"].tolist(),
                      "text":hov_sL,"hoverinfo":"text","showlegend":True})
    traces_id.append({**base,"type":"scatter","name":f"{CLASS_EMOJI[cls]} {cls}",
                      "x":sub["log_bxe"].tolist(),"y":sub["log_s"].tolist(),
                      "text":hov_id,"hoverinfo":"text","showlegend":False})

# Also add ITP points to 2D (s vs L)
itp_sL = []
for sex in ["f","m"]:
    sub = itp[itp["sex"]==sex]
    itp_sL.append({
        "type":"scatter","name":f"ITP {'♀' if sex=='f' else '♂'}",
        "x":np.log10(sub["L"]).tolist(),"y":np.log10(sub["s"].clip(lower=1e-3)).tolist(),
        "mode":"markers","marker":{"color":SEX_COLOR[sex],"symbol":"star","size":8,"opacity":0.9,
                                    "line":{"color":"white","width":0.5}},
        "text":[f"<b>{r['group']}</b> ITP {'♀' if sex=='f' else '♂'}<br>L={r['L']:.0f} d  s={r['s']:.2f}" for _,r in sub.iterrows()],
        "hoverinfo":"text","showlegend":True,"legendgroup":"itp2",
    })
traces_sL.extend(itp_sL)

# ── Phylo violins ──────────────────────────────────────────────────────────
coords  = zims[["lx","ly","lz"]].values.astype(np.float64)
genera  = zims["genus"].values
families= zims["family"].values
orders  = zims["order"].values
cls_arr = zims["class"].values
n = len(zims)
i_idx, j_idx = np.triu_indices(n, k=1)
diff = coords[i_idx]-coords[j_idx]
dists_all = np.sqrt((diff**2).sum(axis=1))
same_gen = genera[i_idx]==genera[j_idx]
same_fam = families[i_idx]==families[j_idx]
same_ord = orders[i_idx]==orders[j_idx]
same_cls = cls_arr[i_idx]==cls_arr[j_idx]
level = np.where(same_gen,0,np.where(same_fam,1,np.where(same_ord,2,np.where(same_cls,3,4))))
lvl_labels = ["Same genus","Same family","Same order","Same class","Diff class"]
lvl_colors = ["#ff3366","#ff8800","#ffdd00","#44cc44","#4499ff"]
rng = np.random.default_rng(42)
phylo_traces = []
for lv in range(5):
    mask = level==lv
    d_lv = dists_all[mask]
    if not len(d_lv): continue
    idx = rng.choice(len(d_lv), min(8000,len(d_lv)), replace=False)
    phylo_traces.append({
        "type":"violin","name":lvl_labels[lv],"y":d_lv[idx].tolist(),
        "box":{"visible":True},"meanline":{"visible":True},
        "fillcolor":lvl_colors[lv],"opacity":0.7,"line":{"color":lvl_colors[lv]},
        "points":False,
    })

# Genus dropdown
genus_counts = zims.groupby("genus").size()
genus_list   = sorted(genus_counts[genus_counts>=3].index)
genus_opts   = '<option value="">— All genera —</option>\n'
for g in genus_list:
    sub = zims[zims["genus"]==g]
    cls0 = sub["class"].iloc[0]
    genus_opts += f'<option value="{g}">{CLASS_EMOJI.get(cls0,"")} {g} ({sub["binSpecies"].nunique()} spp.)</option>\n'

# Class checkbox HTML
cls_cb_html = ""
for cls in classes:
    col = CLASS_COLOR[cls]; em = CLASS_EMOJI[cls]
    cls_cb_html += (
        f'<label class="cb-row" style="border-left:3px solid {col};padding-left:6px;">'
        f'<input type="checkbox" class="cls-cb" value="{cls}" checked> {em} {cls}</label>\n'
    )

# JSON
itp_traces_json    = json.dumps(itp_traces)
zims_traces_json   = json.dumps(zims_traces)
traces_sL_json     = json.dumps(traces_sL)
traces_id_json     = json.dumps(traces_id)
phylo_traces_json  = json.dumps(phylo_traces)
class_color_json   = json.dumps(CLASS_COLOR)
itp_cat_map_json   = itp_cat_map_js
orig_idx_json      = json.dumps(list(range(3,14)))  # existing manifold data indices

# ══════════════════════════════════════════════════════════════════════════════
# Build injection
# ══════════════════════════════════════════════════════════════════════════════
inject = f"""
<style>
/* ── layout ── */
#sr-tabbar {{
  position:fixed;top:0;left:0;right:0;height:42px;z-index:10000;
  background:#111827;border-bottom:2px solid #374151;
  display:flex;align-items:center;padding:0 10px;gap:4px;
}}
.sr-tab {{
  padding:8px 16px;cursor:pointer;background:transparent;border:none;
  color:#9ca3af;font-size:13px;font-weight:600;border-radius:6px;
  transition:all .15s;font-family:Arial,sans-serif;
}}
.sr-tab:hover   {{color:#e5e7eb;background:#374151;}}
.sr-tab.active  {{color:#60a5fa;background:#1f2937;}}
/* push manifold div below tab bar */
#{DIV_ID} {{margin-top:42px;height:calc(100vh - 42px) !important;}}
/* overlay views */
.sr-view {{
  display:none;position:fixed;top:42px;left:0;right:270px;bottom:0;
  z-index:9000;background:#111827;
}}
.sr-view.active {{display:flex;gap:0;}}
#sr-view-phylo {{flex-direction:column;}}
#sr-plot-sL, #sr-plot-id {{flex:1;height:100%;}}
#sr-plot-phylo {{flex:1;}}
#sr-phylo-cap {{
  padding:7px 14px;font-size:10px;color:#6b7280;background:#1f2937;
  border-top:1px solid #374151;font-family:Arial,sans-serif;
}}
/* ── unified filter panel ── */
#sr-panel {{
  position:fixed;top:42px;right:0;bottom:0;width:270px;z-index:9500;
  background:#1f2937;border-left:1px solid #374151;
  overflow-y:auto;padding:12px;font-size:12px;color:#e5e7eb;
  font-family:Arial,sans-serif;
}}
#sr-panel h4 {{margin:0 0 6px;font-size:11px;color:#9ca3af;text-transform:uppercase;letter-spacing:.5px;font-weight:700;}}
.sr-section {{margin-bottom:14px;border-bottom:1px solid #374151;padding-bottom:12px;}}
.sr-section-title {{
  font-weight:700;font-size:13px;color:#60a5fa;margin-bottom:8px;
  display:flex;align-items:center;justify-content:space-between;cursor:pointer;
}}
.sr-section-title span.toggle {{font-size:10px;color:#6b7280;}}
.sr-body {{transition:all .2s;}}
.cb-row {{display:flex;align-items:center;gap:5px;padding:2px 4px;margin-bottom:2px;border-radius:4px;cursor:pointer;}}
.cb-row:hover {{background:#374151;}}
.cb-row input {{margin:0;cursor:pointer;}}
.row2 {{display:grid;grid-template-columns:1fr 1fr;gap:3px;}}
select {{width:100%;padding:4px;background:#374151;border:1px solid #4b5563;
         border-radius:5px;color:#e5e7eb;font-size:11px;margin-bottom:6px;}}
input[type=range] {{width:100%;accent-color:#60a5fa;}}
.rms-val {{color:#fbbf24;font-size:10px;float:right;}}
.stats {{margin-top:2px;font-size:10px;color:#6b7280;line-height:1.7;}}
.stats b {{color:#9ca3af;}}
.good {{color:#34d399;}} .warn {{color:#fbbf24;}} .bad {{color:#f87171;}}
</style>

<!-- Tab bar -->
<div id="sr-tabbar">
  <button class="sr-tab active" onclick="srTab('3d')">🌐 3D Manifold</button>
  <button class="sr-tab"        onclick="srTab('2d')">📊 2D Analysis</button>
  <button class="sr-tab"        onclick="srTab('phylo')">🌳 Phylo</button>
</div>

<!-- 2D overlay -->
<div id="sr-view-2d" class="sr-view">
  <div id="sr-plot-sL"></div>
  <div id="sr-plot-id" style="border-left:1px solid #374151;"></div>
</div>

<!-- Phylo overlay -->
<div id="sr-view-phylo" class="sr-view">
  <div id="sr-plot-phylo"></div>
  <div id="sr-phylo-cap">
    Euclidean distance in log(ρ<sub>β</sub>,ρ<sub>η</sub>,ρ<sub>ε</sub>) space vs taxonomic proximity (GBIF taxonomy).
    Sampled ≤ 8 000 pairs per level from 1.33 M total pairs across 1,633 ZIMS curves.
  </div>
</div>

<!-- Unified filter panel -->
<div id="sr-panel">

  <!-- ITP section -->
  <div class="sr-section">
    <div class="sr-section-title" onclick="toggleSection('itp-body')">
      🐭 NIA ITP <span class="toggle" id="itp-tog">▾</span>
    </div>
    <div class="sr-body" id="itp-body">
      <h4>Intervention</h4>
      <select id="itp-grp">{itp_opts_html}</select>
      <h4>Sex</h4>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="itp-f" checked> <span style="color:#e8634a;">♀ Female</span></label>
        <label class="cb-row"><input type="checkbox" id="itp-m" checked> <span style="color:#2c7bb6;">♂ Male</span></label>
      </div>
      <div class="stats">● removal &nbsp; ◆ production<br>Larger/outlined = Control</div>
    </div>
  </div>

  <!-- ZIMS section -->
  <div class="sr-section">
    <div class="sr-section-title" onclick="toggleSection('zims-body')">
      🦁 ZIMS Zoo Animals <span class="toggle" id="zims-tog">▾</span>
    </div>
    <div class="sr-body" id="zims-body">
      <h4>Class</h4>
      {cls_cb_html}
      <h4 style="margin-top:6px;">Sex</h4>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="zims-f" checked> ♀ Female</label>
        <label class="cb-row"><input type="checkbox" id="zims-m" checked> ♂ Male</label>
      </div>
      <h4 style="margin-top:6px;">Arm</h4>
      <div class="row2">
        <label class="cb-row"><input type="checkbox" id="zims-rem" checked> removal</label>
        <label class="cb-row"><input type="checkbox" id="zims-pro" checked> production</label>
      </div>
      <h4 style="margin-top:6px;">Max RMS <span class="rms-val" id="rms-disp">all</span></h4>
      <input type="range" id="rms-slider" min="0" max="160" value="160" step="1"
             oninput="updateRmsLabel(this.value);applyFilter()">
      <h4 style="margin-top:6px;">Genus</h4>
      <select id="genus-sel" onchange="applyFilter()">{genus_opts}</select>
    </div>
  </div>

  <!-- Legend / notes -->
  <div class="stats">
    <b>Manifold note:</b><br>
    Surface computed at κ=0.5<br>
    ZIMS fits: κ free<br>
    ITP fits:  κ free<br><br>
    <b>ZIMS summary</b><br>
    Curves: <b>{len(zims)}</b><br>
    Median RMS: <b class="good">{zims['rms'].median():.4f}</b><br>
    ndims=5: <b class="warn">{int((zims['ndims']==5).sum())}</b> ({int((zims['ndims']==5).sum())/len(zims)*100:.1f}%)<br>
    Poor fits: <b class="bad">{int((zims['rms']>0.10).sum())}</b><br><br>
    ● ♀ circle &nbsp; ■ ♂ square
  </div>
</div>

<script>
(function(){{
var DIV_ID       = "{DIV_ID}";
var gd           = document.getElementById(DIV_ID);
var itpTraces    = {itp_traces_json};
var zimsTraces   = {zims_traces_json};
var tracesSL     = {traces_sL_json};
var tracesId     = {traces_id_json};
var phyloTr      = {phylo_traces_json};
var catMap       = {itp_cat_map_json};
var clsColor     = {class_color_json};
var ORIG_IDX     = {orig_idx_json};
var itpStart, zimsStart;
var tab2done=false, tabPhyloDone=false;

// ── section collapse ────────────────────────────────────────────────────────
window.toggleSection = function(id) {{
  var el=document.getElementById(id);
  var tog=document.getElementById(id.replace('body','tog'));
  var hidden = el.style.display==='none';
  el.style.display = hidden?'block':'none';
  if(tog) tog.textContent = hidden?'▾':'▸';
}};

// ── tab switching ────────────────────────────────────────────────────────────
window.srTab = function(tab) {{
  document.querySelectorAll('.sr-tab').forEach(function(b,i){{
    b.classList.toggle('active', ['3d','2d','phylo'][i]===tab);
  }});
  ['2d','phylo'].forEach(function(t){{
    document.getElementById('sr-view-'+t).classList.toggle('active', t===tab);
  }});
  if(tab==='2d'   && !tab2done)    {{ init2D();      tab2done=true;     }}
  if(tab==='phylo'&& !tabPhyloDone){{ initPhylo();   tabPhyloDone=true; }}
  if(tab==='3d')   Plotly.Plots.resize(DIV_ID);
}};

// ── dark layout helper ───────────────────────────────────────────────────────
function darkLayout(title,xlab,ylab){{
  return {{
    paper_bgcolor:'#111827',plot_bgcolor:'#1f2937',
    font:{{color:'#e5e7eb',family:'Arial,sans-serif'}},
    xaxis:{{title:xlab,gridcolor:'#374151',zerolinecolor:'#4b5563',color:'#9ca3af'}},
    yaxis:{{title:ylab,gridcolor:'#374151',zerolinecolor:'#4b5563',color:'#9ca3af'}},
    legend:{{bgcolor:'rgba(31,41,55,0.9)',bordercolor:'#4b5563',borderwidth:1,font:{{color:'#e5e7eb',size:11}}}},
    title:{{text:title,font:{{color:'#93c5fd',size:14}}}},
    margin:{{l:65,r:20,t:48,b:60}},
  }};
}}

function init2D(){{
  Plotly.newPlot('sr-plot-sL', tracesSL,
    darkLayout('Sharpness s vs Mean Lifespan L (ZIMS + ITP ★)','log₁₀ L (days)','log₁₀ s'),
    {{responsive:true}});
  Plotly.newPlot('sr-plot-id', tracesId,
    darkLayout('Identifiable Parameter Space','log₁₀(β·Xc/ε)','log₁₀ s'),
    {{responsive:true}});
}}

function initPhylo(){{
  Plotly.newPlot('sr-plot-phylo', phyloTr, {{
    paper_bgcolor:'#111827',plot_bgcolor:'#1f2937',
    font:{{color:'#e5e7eb',family:'Arial,sans-serif'}},
    yaxis:{{title:'Euclidean distance in log ρ-space',gridcolor:'#374151',color:'#9ca3af'}},
    xaxis:{{color:'#9ca3af'}},
    title:{{text:'Parameter-Space Distance vs Taxonomic Proximity',font:{{color:'#93c5fd',size:14}}}},
    violinmode:'overlay',showlegend:true,
    legend:{{bgcolor:'rgba(31,41,55,0.9)',bordercolor:'#4b5563',borderwidth:1,font:{{color:'#e5e7eb',size:11}}}},
    margin:{{l:70,r:20,t:50,b:40}},
  }},{{responsive:true}});
}}

// ── add traces to 3D manifold ────────────────────────────────────────────────
function tryAdd(){{
  if(!gd._fullLayout){{ setTimeout(tryAdd,200); return; }}
  itpStart = gd.data.length;
  Plotly.addTraces(DIV_ID, itpTraces).then(function(){{
    zimsStart = gd.data.length;
    Plotly.addTraces(DIV_ID, zimsTraces).then(function(){{
      // extend show/hide buttons
      var allNew=[];
      for(var i=itpStart;i<gd.data.length;i++) allNew.push(i);
      var ext = ORIG_IDX.concat(allNew);
      Plotly.relayout(DIV_ID,{{
        'updatemenus[0].buttons[0].args':[{{visible:true}},      ext],
        'updatemenus[0].buttons[1].args':[{{visible:'legendonly'}}, ext],
      }});
      applyFilter();
    }});
  }});
}}
tryAdd();

// ── filter helpers ────────────────────────────────────────────────────────────
function itpGroupsFor(val){{
  if(val==='all') return null;
  if(val.startsWith('cat:')) return catMap[val.slice(4)]||[];
  return [val];
}}

window.updateRmsLabel = function(v){{
  document.getElementById('rms-disp').textContent = v>=160 ? 'all' : (v/1000).toFixed(3);
}};

window.applyFilter = function(){{
  if(itpStart===undefined||zimsStart===undefined) return;

  // ITP filter
  var itpGrp  = document.getElementById('itp-grp').value;
  var itpShowF = document.getElementById('itp-f').checked;
  var itpShowM = document.getElementById('itp-m').checked;
  var allowed  = itpGroupsFor(itpGrp);
  var itpVis   = itpTraces.map(function(t){{
    var m=t.meta;
    var sexOk=(m.sex==='f'&&itpShowF)||(m.sex==='m'&&itpShowM);
    var grpOk=(allowed===null)||(allowed.indexOf(m.group)!==-1);
    return sexOk && grpOk;
  }});
  Plotly.restyle(DIV_ID,{{visible:itpVis}}, itpTraces.map(function(_,i){{return itpStart+i;}}));

  // ZIMS filter
  var showCls={{}};
  document.querySelectorAll('.cls-cb').forEach(function(cb){{showCls[cb.value]=cb.checked;}});
  var zShowF = document.getElementById('zims-f').checked;
  var zShowM = document.getElementById('zims-m').checked;
  var zRem   = document.getElementById('zims-rem').checked;
  var zPro   = document.getElementById('zims-pro').checked;
  var rmsMax = document.getElementById('rms-slider').value/1000;
  var genus  = document.getElementById('genus-sel').value;
  var rmsAll = rmsMax>=0.159;

  var zVis=[], zSize=[], zIdx=[];
  zimsTraces.forEach(function(t,i){{
    var m=t.meta;
    if(m.sex==='lbl'){{
      zVis.push(showCls[m.class]?true:false);
      zSize.push([20]); zIdx.push(zimsStart+i); return;
    }}
    var clsOk=(showCls[m.class]);
    var sexOk=(m.sex==='f'&&zShowF)||(m.sex==='m'&&zShowM);
    var armOk=(m.arm==='removal'&&zRem)||(m.arm==='production'&&zPro);
    zVis.push(clsOk&&sexOk&&armOk);
    var sizes=m.rms.map(function(r,j){{
      var rOk=rmsAll||r<=rmsMax;
      var gOk=!genus||m.genus[j]===genus;
      return (rOk&&gOk)?3:0;
    }});
    zSize.push(sizes); zIdx.push(zimsStart+i);
  }});
  Plotly.restyle(DIV_ID,{{visible:zVis}}, zIdx);
  Plotly.restyle(DIV_ID,{{'marker.size':zSize}}, zIdx);
}};

// wire controls
document.querySelectorAll('.cls-cb,#itp-f,#itp-m,#itp-grp,#zims-f,#zims-m,#zims-rem,#zims-pro,#rms-slider').forEach(function(el){{
  el.addEventListener('change', applyFilter);
}});
}})();
</script>
"""

with open(HTML_IN) as f:
    html = f.read()
html = html.replace("</body>", inject + "\n</body>", 1)
with open(HTML_OUT, "w") as f:
    f.write(html)

import os
print(f"Done — {len(itp_traces)} ITP + {len(zims_traces)} ZIMS traces")
print(f"Saved: {HTML_OUT}  ({os.path.getsize(HTML_OUT)/1e6:.1f} MB)")
