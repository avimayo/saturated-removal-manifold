"""
Injects ITP FP-fit points with intervention dropdown + sex toggle
into the saturated-removal manifold HTML.
"""
import pandas as pd, numpy as np, json

ITP_CSV = (
    "/Users/Avimayo/Library/CloudStorage/GoogleDrive-mayoavi@gmail.com"
    "/.shortcut-targets-by-id/1JUcTeGwLojxa9QVWv-SdPTq9kAz4r-xR"
    "/SysMedBook2021/Chapter 6_Aging/ITP survival analysis/itp_fp_fits.csv"
)
HTML_IN  = "/Users/Avimayo/Documents/saturated-removal-manifold/index.html"
HTML_OUT = "/Users/Avimayo/Documents/saturated-removal-manifold/itp_overlay.html"
DIV_ID   = "a1d89836-0daf-4582-ae54-c7c4aca77de6"

df = pd.read_csv(ITP_CSV)
df["lx"] = np.log10(df["rho_beta"])
df["ly"] = np.log10(df["rho_eta"])
df["lz"] = np.log10(df["rho_eps"])

# ── colours & symbols ────────────────────────────────────────────────────────
SEX_COLOR  = {"f": "#e8634a", "m": "#2c7bb6"}
ARM_SYMBOL = {"removal": "circle", "production": "diamond"}

# ── build one trace per (group, sex) ────────────────────────────────────────
groups_order = sorted(df["group"].unique())
traces = []
for grp in groups_order:
    for sex in ["f", "m"]:
        sub = df[(df["group"] == grp) & (df["sex"] == sex)]
        if sub.empty:
            continue
        is_ctrl = (grp == "Control")
        color  = SEX_COLOR[sex]
        symbols = [ARM_SYMBOL[a] for a in sub["arm"]]
        # build hover
        hover = [
            f"<b>{r['group']}</b> | {r['cohort']} | {'♀' if r['sex']=='f' else '♂'}<br>"
            f"L = {r['L']:.0f} d &nbsp; s = {r['s']:.2f} &nbsp; rms = {r['rms']:.3f}<br>"
            f"arm = {r['arm']}<br>"
            f"ρ<sub>β</sub>={r['rho_beta']:.3f}  ρ<sub>η</sub>={r['rho_eta']:.3f}  ρ<sub>ε</sub>={r['rho_eps']:.3f}"
            for _, r in sub.iterrows()
        ]
        trace = {
            "type": "scatter3d",
            "name": f"{grp} ({'♀' if sex=='f' else '♂'})",
            "x": sub["lx"].tolist(),
            "y": sub["ly"].tolist(),
            "z": sub["lz"].tolist(),
            "text": hover,
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {
                "color": color,
                "symbol": symbols,
                "size": 7 if is_ctrl else 5,
                "opacity": 0.95,
                "line": {"color": "black" if is_ctrl else color,
                         "width": 1.5 if is_ctrl else 0},
            },
            "legendgroup": grp,
            "showlegend": True,
            "meta": {"group": grp, "sex": sex},
        }
        traces.append(trace)

traces_json = json.dumps(traces)
groups_json = json.dumps(groups_order)

# ── HTML/CSS/JS injection ────────────────────────────────────────────────────
# Dropdown option groups: map intervention → category
CATEGORIES = {
    "Rapamycin":    ["Rapa_hi","Rapa_lo","Rapa_mid",
                     "Rapa_hi_continuous","Rapa_hi_cycle","Rapa_hi_start_stop"],
    "17α-E2":       ["17aE2","17aE2_hi","17aE2_16m","17aE2_20m"],
    "Metabolic":    ["Met","MetRapa","ACA","Cana"],
    "Antioxidants": ["MitoQ","GGA","MB","bGPA","NR","Min","CC"],
    "Other":        ["HBX","INT-767","DMAG","MIF098","Prot","UDCA"],
    "Control":      ["Control"],
}
# build <optgroup> HTML — with a "select whole category" option per group
optgroups_html = '<option value="all">— All interventions —</option>\n'
for cat, members in CATEGORIES.items():
    optgroups_html += f'<optgroup label="── {cat} ──">\n'
    optgroups_html += f'  <option value="cat:{cat}">★ All {cat}</option>\n'
    for g in members:
        optgroups_html += f'  <option value="{g}">&nbsp;&nbsp;{g}</option>\n'
    optgroups_html += '</optgroup>\n'

# build a JS map from category name → list of group strings (for the filter)
cat_map_js = json.dumps({cat: members for cat, members in CATEGORIES.items()})

inject = f"""
<!-- ITP Control Panel -->
<div id="itp-panel" style="
  position:fixed; top:14px; right:14px; z-index:9999;
  background:rgba(255,255,255,0.95); padding:12px 14px;
  border-radius:10px; box-shadow:0 3px 14px rgba(0,0,0,0.25);
  font-family:Arial,sans-serif; font-size:13px; min-width:195px;">
  <div style="font-weight:bold; font-size:14px; margin-bottom:8px;">
    🐭 NIA ITP overlay
  </div>

  <label style="display:block; margin-bottom:4px; font-weight:600; color:#444;">
    Intervention
  </label>
  <select id="itp-grp" style="width:100%; padding:4px; border-radius:5px;
    border:1px solid #ccc; font-size:12px; margin-bottom:10px;">
    {optgroups_html}
  </select>

  <div style="font-weight:600; color:#444; margin-bottom:5px;">Sex</div>
  <label style="margin-right:12px; cursor:pointer;">
    <input type="checkbox" id="itp-f" checked>
    <span style="color:#e8634a;">♀ Female</span>
  </label>
  <label style="cursor:pointer;">
    <input type="checkbox" id="itp-m" checked>
    <span style="color:#2c7bb6;">♂ Male</span>
  </label>

  <div style="margin-top:10px; font-size:10px; color:#888; line-height:1.5;">
    ● removal arm &nbsp; ◆ production arm<br>
    Larger / outlined = control group<br>
    The "Mice" legend entry has 2 dots:<br>
    ♀ upper-right &nbsp; ♂ lower-left
  </div>
</div>

<script>
(function() {{
  var divId  = "{DIV_ID}";
  var gd     = document.getElementById(divId);
  var traces = {traces_json};
  var catMap = {cat_map_js};
  var itpStart;

  // Original species-data indices the manifest buttons already know about
  var ORIG_DATA_IDX = [3,4,5,6,7,8,9,10,11,12,13];

  function tryAdd() {{
    if (!gd._fullLayout) {{ setTimeout(tryAdd, 200); return; }}
    itpStart = gd.data.length;
    Plotly.addTraces(divId, traces).then(function() {{
      // Extend "show all data" / "hide all data" buttons to cover ITP traces too
      var allNew = [];
      for (var i = itpStart; i < gd.data.length; i++) allNew.push(i);
      var showIdx = ORIG_DATA_IDX.concat(allNew);
      var hideIdx = ORIG_DATA_IDX.concat(allNew);
      Plotly.relayout(divId, {{
        'updatemenus[0].buttons[0].args': [{{'visible': true}},        showIdx],
        'updatemenus[0].buttons[1].args': [{{'visible': 'legendonly'}}, hideIdx],
      }});
      applyFilter();
    }});
  }}
  tryAdd();

  function groupsForSelection(val) {{
    if (val === 'all') return null;
    if (val.startsWith('cat:')) {{
      var cat = val.slice(4);
      return catMap[cat] || [];
    }}
    return [val];
  }}

  function applyFilter() {{
    if (itpStart === undefined) return;
    var sel     = document.getElementById('itp-grp').value;
    var showF   = document.getElementById('itp-f').checked;
    var showM   = document.getElementById('itp-m').checked;
    var allowed = groupsForSelection(sel);

    var visItp = traces.map(function(t) {{
      var sexOk = (t.meta.sex === 'f' && showF) || (t.meta.sex === 'm' && showM);
      var grpOk = (allowed === null) || (allowed.indexOf(t.meta.group) !== -1);
      return sexOk && grpOk;
    }});
    var idxItp = traces.map(function(_, i) {{ return itpStart + i; }});
    Plotly.restyle(divId, {{visible: visItp}}, idxItp);
  }}

  document.getElementById('itp-grp').addEventListener('change', applyFilter);
  document.getElementById('itp-f').addEventListener('change', applyFilter);
  document.getElementById('itp-m').addEventListener('change', applyFilter);
}})();
</script>
"""

with open(HTML_IN) as f:
    html = f.read()

html = html.replace("</body>", inject + "\n</body>", 1)

with open(HTML_OUT, "w") as f:
    f.write(html)

n_traces = len(traces)
print(f"Done — {n_traces} ITP traces injected")
print(f"Saved: {HTML_OUT}")
