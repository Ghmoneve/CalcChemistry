"""Generate every figure in the paper from the JSON artefacts in experiments/.

WHY THIS IS A SEPARATE SCRIPT
Every number plotted here already exists on disk. Regenerating the figures must not
require rerunning a DFT or an FCI calculation, so this script reads only JSON and
never imports pyscf. If a figure cannot be drawn from stored data, the honest fix is
to store the data, not to recompute it inside the plotting code.

WHAT EACH FIGURE CARRIES
Readers should be able to check each panel against a table:
  fig_tau_estimates    the four estimates of tau against the diagonal gap
  fig_sector_trend     tau and d as the reference sector grows
  fig_route_slack      slack of the block, Temple and Schur certificates
  fig_similarity_floor how far random diagonal scalings can push the norm
  fig_schur_sensitivity  the validity cliff in the assumed spectral floor
  fig_knobs            cost of the shift and reference-block knobs
  fig_caf2_bend        bending potential, DFT against CCSD(T) and experiment
  fig_barrier_bar      bend barrier across levels
  fig_reaction_bracket reaction energy across eight levels, two treatments
  fig_grid_convergence imaginary mode against DFT integration grid
  fig_coverage         growth per cycle mapped to fractional site coverage
  fig_imaginary        spurious imaginary modes found before reoptimisation

Output: paper/figures/*.pdf, vector, sized for a single column.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'paper' / 'figures'
KCAL = 627.5094740631

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 8.0,
    'axes.labelsize': 8.0,
    'axes.titlesize': 8.5,
    'legend.fontsize': 7.5,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'figure.dpi': 200,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
})

C_MAIN = '#1f4e79'
C_ACC = '#c55a11'
C_GREY = '#7f7f7f'
C_GREEN = '#2e7d32'


def load(rel):
    p = ROOT / rel
    if not p.exists():
        print(f'  MISSING {rel}')
        return None
    return json.loads(p.read_text())


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / f'{name}.pdf'
    fig.savefig(out)
    if os.environ.get('FIG_PNG'):
        fig.savefig(FIG / f'{name}.png', dpi=150)
    plt.close(fig)
    print(f'  wrote {out.relative_to(ROOT)}')


# ---------------------------------------------------------------- figure 1
def fig_tau_estimates():
    d = load('experiments/complement/tau_measure.json')
    if not d:
        return
    names = ['$\\|O\\|_2$\n(spectral)', 'Gershgorin', 'max row\n2-norm', 'Frobenius']
    vals = [d['tau_spectral'], d['tau_gershgorin'],
            d['tau_max_row_2norm'], d['tau_frobenius']]
    gap = d['diagonal_gap_to_next']
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    cols = [C_MAIN] + [C_GREY] * 3
    bars = ax.bar(names, vals, color=cols, width=0.62, zorder=2)
    ax.axhline(gap, color=C_ACC, ls='--', lw=1.2, zorder=3,
               label=f'diagonal gap $\\Delta={gap:.3f}$ Ha')
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.18, f'{v:.3f}',
                ha='center', va='bottom', fontsize=7.5)
    ax.set_ylabel('$\\tau$  (Ha)')
    ax.set_ylim(0, max(vals) * 1.18)
    ax.legend(loc='upper left', frameon=False)
    ax.set_title('only the spectral norm is sharp')
    save(fig, 'fig_tau_estimates')


# ---------------------------------------------------------------- figure 2
def fig_sector_trend():
    d = load('experiments/complement/tau_measure.json')
    if not d or 'sector_trend' not in d:
        return
    rows = d['sector_trend']
    k = [r['k_active'] for r in rows]
    tau = [r['tau'] for r in rows]
    gap = d['diagonal_gap_to_next']
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    ax.plot(k, tau, 'o-', color=C_MAIN, ms=4, lw=1.2, label='$\\tau=\\|O\\|_2$')
    ax.axhline(gap, color=C_ACC, ls='--', lw=1.2,
               label=f'diagonal gap ${gap:.3f}$ Ha')
    ax.set_xlabel('orbitals kept in the reference sector')
    ax.set_ylabel('$\\tau$  (Ha)')
    ax.set_yscale('log')
    ax.legend(frameon=False)
    ax.set_title('the gap does not shrink with the sector')
    save(fig, 'fig_sector_trend')


# ---------------------------------------------------------------- figure 3
def fig_route_slack():
    """Three rungs of the same certificate, each using more information.

    The numbers are recomputed from stored primitives rather than read from
    route_comparison.json, whose temple_true_gap entries were produced with an
    inconsistent residual and yield bounds ABOVE E_FCI, i.e. not bounds at all.
    Each value plotted here is asserted to lie below E_FCI.
    """
    sr = load('experiments/complement/schur_route.json')
    tm = load('experiments/complement/tau_measure.json')
    vc = load('experiments/complement/verified_cholesky.json')
    if not (sr and tm and vc):
        return
    d, efci, b = sr['d'], sr['E_FCI'], sr['beta_row']
    delta = sr['lambda_min_Hr'] - d
    slack_norm = (d - tm['tau_spectral'] - efci) * KCAL
    scalar = d - 2 * b * b / (np.sqrt(delta ** 2 + 4 * b * b) + delta)
    slack_scalar = (scalar - efci) * KCAL
    slack_schur = (vc['schur_bound_certified'] - efci) * KCAL
    for name, v in (('full norm', slack_norm), ('scalar', slack_scalar),
                    ('schur', slack_schur)):
        assert v < 0, f'{name} bound lies above E_FCI and is not a bound'
    labels = ['full norm\n$d-\\|O\\|_2$',
              'one row + gap\n$d-\\frac{2b^2}{\\sqrt{\\delta^2+4b^2}+\\delta}$',
              'exact Schur\n$d-b^{\\top}A^{-1}b$']
    vals = [abs(slack_norm), abs(slack_scalar), abs(slack_schur)]
    fig, ax = plt.subplots(figsize=(6.3, 4.26))
    bars = ax.bar(labels, vals, color=[C_GREY, '#5b8db8', C_MAIN], width=0.62, zorder=2)
    ax.set_yscale('log')
    ax.set_ylabel('slack below $E_{\\mathrm{FCI}}$  (kcal/mol)')
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v * 1.18, f'{v:.1f}',
                ha='center', va='bottom', fontsize=7.5)
    ax.set_ylim(0.3, max(vals) * 5)
    for i in (1, 2):
        ax.text(i, vals[i] * 2.6, f'$\\times${vals[i-1]/vals[i]:.0f} tighter',
                color=C_ACC, fontsize=7.0, ha='center')
    ax.set_title('each rung buys about two orders of magnitude')
    save(fig, 'fig_route_slack')


# ---------------------------------------------------------------- figure 4
def fig_similarity_floor():
    au = load('experiments/complement/sign_scaling_audit.json')
    tm = load('experiments/complement/tau_measure.json')
    if not (au and tm):
        return
    gap = tm['diagonal_gap_to_next']
    labels = ['gap $\\Delta$', 'best over\n20k scalings', 'Gershgorin\n(best scaling)']
    vals = [gap, au['brute_force_min_sigma_max_over_20000_random_D'],
            au['brute_force_min_gershgorin_over_2000_random_scalings']]
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    bars = ax.bar(labels, vals, color=[C_ACC, C_MAIN, C_GREY], width=0.6, zorder=2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.06, f'{v:.3f}',
                ha='center', va='bottom', fontsize=7.5)
    ax.set_ylabel('achieved $\\tau$  (Ha)')
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_title('no diagonal scaling reaches the gap')
    save(fig, 'fig_similarity_floor')


# ---------------------------------------------------------------- figure 5
def fig_schur_sensitivity():
    d = load('experiments/complement/schur_route.json')
    if not d or 'spectral_bound_sensitivity' not in d:
        return
    rows = d['spectral_bound_sensitivity']
    mu = [r['mu'] for r in rows]
    g = [r['g_mu'] for r in rows]
    efci = d['E_FCI']
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    ax.plot(mu, g, 'o-', color=C_MAIN, ms=4, lw=1.2, label='$g(\\mu)=d-b^{\\top}(H_r-\\mu I)^{-1}b$')
    ax.axhline(efci, color=C_ACC, ls='--', lw=1.2, label='$E_{\\mathrm{FCI}}$')
    ax.set_xlabel('assumed floor $\\mu$  (Ha)')
    ax.set_ylabel('evaluated scalar  (Ha)')
    ax.legend(frameon=False, loc='best')
    ax.set_title('the bound is only valid left of the cliff')
    save(fig, 'fig_schur_sensitivity')


# ---------------------------------------------------------------- figure 6
def fig_knobs():
    d = load('experiments/complement/knob_tradeoff.json')
    if not d:
        return
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.40))
    ax = axes[0]
    rows = d['shift_scan']
    ax.plot([r['t'] for r in rows], [r['slack_kcal'] for r in rows],
            'o-', color=C_MAIN, ms=4, lw=1.2)
    ax.set_xlabel('shift $t$  (Ha)')
    ax.set_ylabel('slack  (kcal/mol)')
    ax.set_title('(a) shifting the reference level')
    ax2 = axes[1]
    rows = d['reference_space_scan']
    ax2.plot([r['k'] for r in rows], [r['slack_kcal'] for r in rows],
             's-', color=C_ACC, ms=4, lw=1.2)
    ax2.set_xlabel('reference determinants $k$')
    ax2.set_ylabel('slack  (kcal/mol)')
    ax2.set_title('(b) enlarging the reference block')
    fig.tight_layout()
    save(fig, 'fig_knobs')


# ---------------------------------------------------------------- figure 7
def fig_caf2_bend():
    fm = load('experiments/caf2/flat_mode.json')
    cc = load('experiments/caf2/cct_barrier.json')
    if not fm:
        return
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.50))
    ax = axes[0]
    keys = sorted(k for k in fm if k.startswith('scan_'))
    for i, k in enumerate(keys):
        s = fm[k]
        rows = s['rows']
        # dE_kcal in flat_mode.json is measured from each level's own minimum, so
        # theta=180 sits at +2.2 kcal. CCSD(T) rel_kcal is measured from the linear
        # geometry. Both panels must use E(theta) - E(180) or the contrast inverts.
        e_lin = next((r['dE_kcal'] for r in rows if r['theta'] == 180), 0.0)
        th = [r['theta'] for r in rows]
        de = [r['dE_kcal'] - e_lin for r in rows]
        lab = k.replace('scan_', '').replace('def2-', '').replace('_', ' ').upper()
        ax.plot(th, de, 'o-', ms=3, lw=1.0,
                color=plt.cm.viridis(i / max(1, len(keys) - 1)), label=lab)
    ax.set_xlabel('F--Ca--F angle (deg)')
    ax.set_ylabel('$E(\\theta)-E(180^\\circ)$  (kcal/mol)')
    ax.set_title('(a) DFT: an apparent well')
    # every curve now lies at or below zero, so the strip just under the top axis
    # is empty and the legend cannot sit on data
    ax.set_ylim(top=0.32)
    ax.legend(frameon=False, fontsize=7.0, ncol=3, loc='upper left')
    ax2 = axes[1]
    if cc:
        for j, (k, s) in enumerate(sorted(cc.items())):
            ax2.plot(s['angles'], s['rel_kcal'], 'o-', ms=3, lw=1.1,
                     color=plt.cm.plasma(0.20 + 0.55 * j / max(1, len(cc) - 1)),
                     label='CCSD(T) ' + k.replace('def2-', '').replace('_', ' '))
    ax2.plot([154.3], [0.15], '*', ms=11, color='k', zorder=5, label='experiment')
    ax2.axhline(0, color='k', lw=0.5, ls=':')
    ax2.set_xlabel('F--Ca--F angle (deg)')
    ax2.set_ylabel('$\\Delta E$ vs linear (kcal/mol)')
    ax2.set_title('(b) CCSD(T): essentially flat')
    ax2.legend(frameon=False, fontsize=7.0, loc='upper right')
    fig.tight_layout()
    save(fig, 'fig_caf2_bend')


# ---------------------------------------------------------------- figure 8
def fig_barrier_bar():
    fm = load('experiments/caf2/flat_mode.json')
    cc = load('experiments/caf2/cct_barrier.json')
    if not fm:
        return
    labels, vals, cols = [], [], []
    for k in sorted(x for x in fm if x.startswith('scan_')):
        labels.append('DFT  ' + k.replace('scan_', '').replace('def2-', '').replace('_', '/'))
        vals.append(fm[k]['bend_rise_kcal'])
        cols.append(C_GREY)
    if cc:
        for k in sorted(cc):
            labels.append('CCSD(T)  ' + k.replace('def2-', '').replace('_', ' '))
            vals.append(cc[k]['barrier_kcal'])
            cols.append(C_MAIN)
    labels.append('experiment')
    vals.append(0.15)
    cols.append(C_ACC)
    order = np.argsort(vals)
    labels = [labels[i] for i in order]
    cols = [cols[i] for i in order]
    vals = [vals[i] for i in order]
    fig, ax = plt.subplots(figsize=(6.3, 4.73))
    y = np.arange(len(vals))
    ax.barh(y, vals, color=cols, height=0.66, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.0)
    ax.set_xlabel('bend barrier (kcal/mol)')
    ax.axvline(0.15, color=C_ACC, ls=':', lw=1.0)
    for yi, v in zip(y, vals):
        ax.text(v + 0.05, yi, f'{v:.2f}', va='center', fontsize=7.0)
    ax.set_xlim(0, max(vals) * 1.22)
    ax.set_title('DFT overestimates the barrier by up to 20x')
    save(fig, 'fig_barrier_bar')


# ---------------------------------------------------------------- figure 9
def fig_reaction_bracket():
    r = load('experiments/caf2/reaction3.json')
    if not r:
        return
    levels = sorted(r['single_points'])
    bent = [r['single_points'][k]['dE_bent_kcal'] for k in levels]
    lin = [r['single_points'][k]['dE_linear_kcal'] for k in levels]
    s = r['summary']
    fig, ax = plt.subplots(figsize=(6.3, 4.55))
    x = np.arange(len(levels))
    ax.axhspan(s['bracket_low'], s['bracket_high'], color=C_GREY, alpha=0.16, zorder=0,
               label='bracket from all levels')
    ax.plot(x, bent, 'o', color=C_MAIN, ms=5, label='bent CaF$_2$')
    ax.plot(x, lin, 's', color=C_ACC, ms=4, label='linear CaF$_2$')
    ax.set_xticks(x)
    short = [l.replace('def2-', '').replace('b3lyp', 'B3LYP').replace('pbe', 'PBE')
             for l in levels]
    ax.set_xticklabels(short, rotation=40, ha='right', fontsize=7.0)
    ax.set_ylabel('reaction energy  (kcal/mol)')
    ax.legend(frameon=False, fontsize=7.5)
    ax.set_title(f"spread {s['dE_bent_spread']:.2f} kcal/mol across levels")
    save(fig, 'fig_reaction_bracket')


# ---------------------------------------------------------------- figure 10
def fig_grid_convergence():
    fm = load('experiments/caf2/flat_mode.json')
    if not fm or 'grid' not in fm:
        return
    rows = fm['grid']
    lv = [r['level'] for r in rows]
    cm = [r['bend_cm'] for r in rows]
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    ax.plot(lv, cm, 'o-', color=C_MAIN, ms=5, lw=1.2)
    ax.axhline(0, color='k', lw=0.6, ls=':')
    for x, y in zip(lv, cm):
        if x <= 3:
            ax.annotate(f'{y:.0f}', (x, y), textcoords='offset points',
                        xytext=(0, 6), ha='center', fontsize=7.0, color=C_GREY)
    ax.set_xlabel('DFT integration grid level')
    ax.set_ylabel('bend mode  (cm$^{-1}$)')
    ax.set_title('the imaginary mode survives grid refinement')
    save(fig, 'fig_grid_convergence')


# ---------------------------------------------------------------- figure 11
def fig_coverage():
    g = load('experiments/caf2/gpc_coverage.json')
    if not g:
        return
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    x = np.array(g['gpc'])
    y = np.array(g['coverage']) * 100
    ax.plot(x, y, 'o-', color=C_MAIN, ms=6, lw=1.4)
    for xi, yi, ai in zip(x, y, g['area_per_metal']):
        ax.annotate(f'{yi:.2f}%  ({ai:.0f} $\\AA^2$/Ca)', (xi, yi),
                    textcoords='offset points', xytext=(6, -10), fontsize=7.0)
    ax.set_xlabel('growth per cycle  ($\\AA$)')
    ax.set_ylabel('fractional Ca site coverage  (%)')
    ax.set_xlim(0.26, 0.44)
    ax.set_title('0.3--0.4 $\\AA$/cycle means sub-monolayer')
    save(fig, 'fig_coverage')


# ---------------------------------------------------------------- figure 12
def fig_imaginary():
    f = load('experiments/caf2/frequencies.json')
    if not f:
        return
    names, counts, mags = [], [], []
    for k, v in f.items():
        names.append(k.replace('_', ' '))
        counts.append(v['imaginary_count'])
        mags.append(abs(min(v['imaginary_cm'])) if v['imaginary_cm'] else 0.0)
    fig, ax = plt.subplots(figsize=(6.3, 4.39))
    x = np.arange(len(names))
    bars = ax.bar(x, mags, color=[C_ACC if c else C_GREEN for c in counts],
                  width=0.6, zorder=2)
    for b, c, m in zip(bars, counts, mags):
        ax.text(b.get_x() + b.get_width() / 2, m + 4,
                f'{c} imaginary' if c else 'minimum',
                ha='center', va='bottom', fontsize=7.0)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha='right', fontsize=7.5)
    ax.set_ylabel('largest imaginary mode  (cm$^{-1}$)')
    ax.set_title('three of four geometries were saddle points')
    save(fig, 'fig_imaginary')


def main():
    print('generating figures from stored JSON')
    for fn in (fig_tau_estimates, fig_sector_trend, fig_route_slack,
               fig_similarity_floor, fig_schur_sensitivity, fig_knobs,
               fig_caf2_bend, fig_barrier_bar, fig_reaction_bracket,
               fig_grid_convergence, fig_coverage, fig_imaginary):
        try:
            fn()
        except Exception as exc:                     # noqa: BLE001
            print(f'  FAILED {fn.__name__}: {type(exc).__name__}: {exc}')


if __name__ == '__main__':
    main()
