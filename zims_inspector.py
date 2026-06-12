"""
ZIMS SR fit inspector — Streamlit app.
Deployable version: uses pre-computed KM cache, no raw data needed.
"""
import sys, math, json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

sys.path.insert(0, 'scripts')

DATA_DIR  = 'data'
FITS_CSV  = 'results/zims_all.csv'
FLAGS_CSV = f'{DATA_DIR}/zims_quality_flags.csv'
KS_CSV    = f'{DATA_DIR}/zims_ks.csv'
KM_CACHE  = f'{DATA_DIR}/zims_km_cache.parquet'

CLASSES = ['Amphibia','Aves','Chondrichthyes','Mammalia','Reptilia']
TAU     = np.linspace(0, 3.0, 600)

CLASS_COLOR = {
    'Mammalia':      '#e07b39',
    'Aves':          '#5b8fce',
    'Reptilia':      '#5bce7a',
    'Amphibia':      '#ce5b9a',
    'Chondrichthyes':'#cece5b',
}

# ── Data loading ──────────────────────────────────────────────────────────────
@st.cache_data
def load_fits():
    fits  = pd.read_csv(FITS_CSV)
    flags = pd.read_csv(FLAGS_CSV)
    ks    = pd.read_csv(KS_CSV)[['class','binSpecies','sex','ks_D','neg_log10_p']]
    fits  = fits.merge(flags[['binSpecies','sex','n_flags','flags']], on=['binSpecies','sex'], how='left')
    fits  = fits.merge(ks, on=['class','binSpecies','sex'], how='left')
    return fits

@st.cache_data
def load_km_cache():
    df = pd.read_parquet(KM_CACHE)
    km = {}
    for _, row in df.iterrows():
        key = (row['class'], row['binSpecies'], row['sex'])
        km[key] = {
            'km_all_tau': json.loads(row['km_all_tau']),
            'km_all_s':   json.loads(row['km_all_s']),
            'km_tau':     json.loads(row['km_tau']),
            'km_s':       json.loads(row['km_s']),
            'na_tau':     json.loads(row['na_tau']),
            'na_H':       json.loads(row['na_H']),
        }
    return km

# ── SR curve ──────────────────────────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def compute_sr(rho_eta, rho_beta, rho_eps, kappa, L, ndims, ext_haz):
    from fpt_full_model import fpt_full_survival
    s_sr, _, _, _ = fpt_full_survival(rho_eta, rho_beta, rho_eps, kappa, TAU)
    if ndims == 5 and pd.notna(ext_haz):
        h    = math.exp(-ext_haz)
        s_sr = s_sr * np.exp(-h * L * TAU)
        s_sr = s_sr / max(s_sr[0], 1e-9)
    sr_chf = -np.log(np.maximum(s_sr, 1e-10))
    return TAU.tolist(), s_sr.tolist(), sr_chf.tolist()

# ── Main plot ─────────────────────────────────────────────────────────────────
def build_figure(selected_rows, km_cache):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Survival  S(τ)', 'Cumulative Hazard  H(τ) = −log S'),
        vertical_spacing=0.12, row_heights=[0.55, 0.45],
    )
    for i, row in enumerate(selected_rows):
        color = CLASS_COLOR.get(row['class'], '#888')
        sp, sx = row['binSpecies'], row['sex']
        lbl = f"{sp} {'♀' if sx=='f' else '♂'}"
        key = (row['class'], sp, sx)

        km = km_cache.get(key, {})
        tau_fit, s_fit, sr_chf = compute_sr(
            row['rho_eta'], row['rho_beta'], row['rho_eps'], row['kappa'],
            row['L'], int(row['ndims']), row.get('external_hazard', float('nan'))
        )

        r, g, b = int(color[1:3],16), int(color[3:5],16), int(color[5:7],16)
        faint = f'rgba({r},{g},{b},0.3)'

        if km:
            fig.add_trace(go.Scatter(
                x=km['km_all_tau'], y=km['km_all_s'], mode='lines',
                line=dict(color=faint, width=1, dash='dot'),
                name='KM (all)', legendgroup=lbl, showlegend=False,
            ), row=1, col=1)
            fig.add_trace(go.Scatter(
                x=km['km_tau'], y=km['km_s'], mode='lines',
                line=dict(color=color, width=2),
                name=lbl, legendgroup=lbl, showlegend=True,
            ), row=1, col=1)
            fig.add_trace(go.Scatter(
                x=km['na_tau'], y=km['na_H'], mode='lines',
                line=dict(color=color, width=2),
                name=lbl, legendgroup=lbl, showlegend=False,
            ), row=2, col=1)

        fig.add_trace(go.Scatter(
            x=tau_fit, y=s_fit, mode='lines',
            line=dict(color=color, width=2, dash='dash'),
            name=f'{lbl} (fit)', legendgroup=lbl, showlegend=True,
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=tau_fit, y=sr_chf, mode='lines',
            line=dict(color=color, width=2, dash='dash'),
            name=f'{lbl} (fit)', legendgroup=lbl, showlegend=False,
        ), row=2, col=1)

    fig.update_xaxes(title_text='τ = t / L', gridcolor='#eee', zeroline=False,
                     title_font=dict(color='#333'), tickfont=dict(color='#333'))
    fig.update_yaxes(gridcolor='#eee', zeroline=False, autorange=True,
                     title_font=dict(color='#333'), tickfont=dict(color='#333'))
    fig.add_vline(x=1.0, line=dict(color='#aaa', width=1, dash='dot'),
                  annotation_text='τ=1', annotation_position='top right',
                  annotation_font=dict(color='#555'))
    fig.update_layout(
        height=700, plot_bgcolor='#fff', paper_bgcolor='#fff',
        font=dict(color='#333'),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0,
                    font=dict(color='#333')),
        margin=dict(l=60, r=20, t=80, b=40),
    )
    return fig

# ── Overview plots ────────────────────────────────────────────────────────────
def build_ks_scatter(df):
    import math as _m
    bonf = -_m.log10(0.05 / len(df))
    traces = []
    for cls in CLASSES:
        sub = df[df['class'] == cls].dropna(subset=['ks_D','neg_log10_p'])
        if sub.empty:
            continue
        traces.append(go.Scatter(
            x=sub['ks_D'], y=sub['neg_log10_p'], mode='markers',
            marker=dict(
                color=CLASS_COLOR.get(cls, '#888'),
                size=(sub['n_dead'].clip(0,2000) / 2000 * 14 + 4).tolist(),
                opacity=0.65, line=dict(width=0.5, color='white'),
            ),
            name=cls,
            text=(sub['binSpecies'] + ' ' + sub['sex'].map({'f':'♀','m':'♂'})
                  + '<br>n_dead=' + sub['n_dead'].astype(str)
                  + '  rms=' + sub['rms'].round(4).astype(str)),
            hovertemplate='%{text}<br>D=%{x:.3f}  −log₁₀p=%{y:.1f}<extra></extra>',
        ))
    fig = go.Figure(traces)
    fig.add_hline(y=bonf, line=dict(color='#d62728', width=1.5, dash='dash'),
                  annotation_text=f'Bonferroni  p=0.05/{len(df)}',
                  annotation_position='top right')
    fig.update_layout(
        title=dict(text='KS effect size D vs −log₁₀(p)  [marker size ∝ n_dead]',
                   font=dict(color='#333')),
        xaxis=dict(title='KS D', title_font=dict(color='#333'), tickfont=dict(color='#333'),
                   gridcolor='#eee', zeroline=False),
        yaxis=dict(title='−log₁₀(p)', title_font=dict(color='#333'), tickfont=dict(color='#333'),
                   gridcolor='#eee'),
        height=450, plot_bgcolor='#fff', paper_bgcolor='#fff',
        font=dict(color='#333'),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, font=dict(color='#333')),
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig

def build_rms_box(df):
    fig = go.Figure()
    for cls in CLASSES:
        sub = df[df['class'] == cls]
        if sub.empty:
            continue
        fig.add_trace(go.Box(
            y=sub['rms'], name=cls, marker_color=CLASS_COLOR.get(cls,'#888'),
            boxpoints='outliers', jitter=0.3,
        ))
    fig.update_layout(
        title=dict(text='RMS by class', font=dict(color='#333')),
        yaxis=dict(title='RMS', title_font=dict(color='#333'), tickfont=dict(color='#333'),
                   gridcolor='#eee'),
        xaxis=dict(tickfont=dict(color='#333')),
        height=380, plot_bgcolor='#fff', paper_bgcolor='#fff',
        font=dict(color='#333'),
        showlegend=False, margin=dict(l=60, r=20, t=50, b=40),
    )
    return fig

def build_flag_bar(df):
    types = [('ok','Clean'),('low_n_dead','Low n'),('poor_rms','Poor RMS'),
             ('infant_spike','Infant spike'),('bimodal_BC','Bimodal'),
             ('high_makeham','High Makeham'),('hazard_bathtub','Bathtub')]
    colors = ['#4caf50','#90caf9','#ef5350','#ffb74d','#ce93d8','#80cbc4','#f48fb1']
    labels, vals = [], []
    for (key, lbl), _ in zip(types, colors):
        labels.append(lbl)
        vals.append((df['flags']=='ok').sum() if key=='ok'
                    else df['flags'].str.contains(key,na=False).sum())
    fig = go.Figure(go.Bar(x=labels, y=vals, marker_color=colors,
                           text=vals, textposition='outside',
                           textfont=dict(color='#333')))
    fig.update_layout(
        title=dict(text=f'Quality flags  (N={len(df)} curves)', font=dict(color='#333')),
        yaxis=dict(title='Count', title_font=dict(color='#333'), tickfont=dict(color='#333'),
                   gridcolor='#eee', range=[0, max(vals)*1.15]),
        xaxis=dict(tickfont=dict(color='#333')),
        height=350, plot_bgcolor='#fff', paper_bgcolor='#fff',
        font=dict(color='#333'),
        showlegend=False,
        margin=dict(l=60, r=20, t=50, b=50),
    )
    return fig

# ── App ───────────────────────────────────────────────────────────────────────
st.set_page_config(page_title='ZIMS SR Inspector', layout='wide')
st.title('ZIMS SR fit inspector')

fits    = load_fits()
km_cache = load_km_cache()

with st.sidebar:
    st.header('Filters')
    classes = st.multiselect('Class', CLASSES, default=CLASSES)
    sexes   = st.multiselect('Sex', ['f','m'], default=['f','m'],
                             format_func=lambda x: '♀ Female' if x=='f' else '♂ Male')
    arms    = st.multiselect('Arm', ['removal','production'], default=['removal','production'])
    quality = st.radio('Quality', ['All', 'Clean only', 'Flagged only'], index=0)
    n_min, n_max = int(fits['n_dead'].min()), int(fits['n_dead'].max())
    n_range = st.slider('n_dead', n_min, n_max, (n_min, n_max))
    rms_max = st.slider('Max RMS', 0.005, float(fits['rms'].max()), float(fits['rms'].max()), step=0.005)
    search  = st.text_input('Species search', '')

mask = (
    fits['class'].isin(classes) & fits['sex'].isin(sexes) &
    fits['arm'].isin(arms) & fits['n_dead'].between(*n_range) &
    (fits['rms'] <= rms_max)
)
if quality == 'Clean only':
    mask &= (fits['n_flags'] == 0)
elif quality == 'Flagged only':
    mask &= (fits['n_flags'] > 0)
if search:
    mask &= fits['binSpecies'].str.contains(search, case=False, na=False)

filtered = fits[mask].copy()
filtered['label'] = (
    filtered['flags'].apply(lambda f: '⚠ ' if (f != 'ok' and pd.notna(f)) else '')
    + filtered['binSpecies'] + '  '
    + filtered['sex'].map({'f':'♀','m':'♂'}) + '  ['
    + filtered['class'].str[:3] + ']'
    + '  rms=' + filtered['rms'].round(4).astype(str)
    + '  n=' + filtered['n_dead'].astype(str)
)

st.caption(f'{len(filtered)} curves  |  {(filtered["n_flags"]==0).sum()} clean, '
           f'{(filtered["n_flags"]>0).sum()} flagged')

tab_inspect, tab_overview = st.tabs(['🔬 Inspect curves', '📊 Dataset overview'])

# ── Inspect ───────────────────────────────────────────────────────────────────
with tab_inspect:
    col1, col2 = st.columns([2, 1])
    with col1:
        options = filtered['label'].tolist()
        if 'chosen' not in st.session_state or not st.session_state.chosen:
            best = filtered.sort_values('rms')['label'].iloc[0] if len(filtered) else None
            st.session_state.chosen = [best] if best else []
        b1, b2 = st.columns([4, 1])
        with b2:
            if st.button('🎲 Random', use_container_width=True):
                st.session_state.chosen = [filtered.sample(1)['label'].iloc[0]]
        with b1:
            chosen = st.multiselect('Select curves (up to 5)', options,
                                    default=[o for o in st.session_state.chosen if o in options],
                                    max_selections=5, key='chosen_select')
        st.session_state.chosen = chosen

    if chosen:
        sel_rows = filtered[filtered['label'].isin(chosen)].to_dict('records')
        with col2:
            st.subheader('Parameters')
            disp_cols = ['binSpecies','sex','class','n_dead','rms','ks_D','neg_log10_p',
                         'ndims','arm','s','beta_xc_eps','kappa','L','n_flags','flags']
            disp_cols = [c for c in disp_cols if c in filtered.columns]
            tbl = (filtered[filtered['label'].isin(chosen)][disp_cols]
                   .set_index('binSpecies').T.astype(str))
            st.dataframe(tbl, width='stretch')
        with st.spinner('Computing SR fits…'):
            fig = build_figure(sel_rows, km_cache)
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info('Select a curve above to inspect it.')
        disp = filtered[['binSpecies','sex','class','n_dead','rms','ks_D',
                         'neg_log10_p','arm','ndims','n_flags','flags']].copy()
        disp['sex'] = disp['sex'].map({'f':'♀','m':'♂'})
        st.dataframe(disp.sort_values('rms').reset_index(drop=True), height=400, width='stretch')

# ── Overview ──────────────────────────────────────────────────────────────────
with tab_overview:
    oc1, oc2 = st.columns(2)
    with oc1:
        st.plotly_chart(build_ks_scatter(filtered), use_container_width=True)
    with oc2:
        st.plotly_chart(build_flag_bar(filtered), use_container_width=True)
    st.plotly_chart(build_rms_box(filtered), use_container_width=True)

    st.subheader('Summary by class')
    summary = (filtered.groupby('class').agg(
        N=('rms','count'),
        clean=('n_flags', lambda x: (x==0).sum()),
        rms_med=('rms','median'),
        rms_p90=('rms', lambda x: x.quantile(0.9)),
        ks_D_med=('ks_D','median'),
        n_dead_med=('n_dead','median'),
    ).reset_index().sort_values('N', ascending=False))
    summary['% clean'] = (summary['clean']/summary['N']*100).round(1).astype(str)+'%'
    st.dataframe(summary.set_index('class').round(4), width='stretch')
