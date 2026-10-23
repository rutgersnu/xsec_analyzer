"""Python reference for the xsec_analyzer nuNCQE implementation.

Uses the selection and weight conventions of NCQE_XS_analysis_v3 (lib/selection.py, lib/extraction.py,
lib/reweight.py; read-only) on exactly the files given, so the C++ chain (annie_ncqe_prep -> univmake ->
ncqe_xsec) can be compared number by number.

  select     event-by-event comparison of the selection with a merged annie_ncqe_prep output
  reference  self-contained extraction and uncertainties (same JSON layout as ncqe_xsec)
  compare    compare two JSON summaries (ncqe_xsec vs reference)

Examples
  python3 reference_python.py select    --mc-dir DIR --list files.txt --prep merged.ncqe.root
  python3 reference_python.py reference --mc-dir DIR --list files.txt --mc-pot 4.41733e20 --out ref.json [--unisims]
  python3 reference_python.py compare   ncqe_xsec.json ref.json
"""
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import uproot

V3 = os.environ.get('NCQE_V3_DIR', os.path.expanduser('~/hep/ANNIE/NCQE_XS_analysis_v3'))
if not os.path.isdir(V3):
    V3 = '/workspace/hep/ANNIE/NCQE_XS_analysis_v3'
sys.path.insert(0, V3)
from lib import selection as SEL          # noqa: E402
from lib import extraction as X           # noqa: E402
from lib import reweight as RW            # noqa: E402

# Data-side constants of the PRD draft (same values as ncqe/ncqe_inputs.conf)
INPUTS = dict(n_obs=772, d_cit=81.29333333333334, d_cit_err=9.93156136320524, d_sky=13.952588005387733,
              d_sky_err=12.868699904856147, pot_data=2.572970591694e20)
LABELS = {'sig': '#nu NCQE (FV)', 'OOFV_nu': '#nu NCQE (OOFV)', 'NCother': 'NC other', 'CC': 'CC',
          'nubar': '#bar{#nu} NCQE', 'OOFV_other': 'OOFV other', 'Ext': 'External'}


def read_list(path):
    with open(path) as f:
        return [l.split()[0] for l in f if l.strip() and not l.startswith('#')]


# ------------------------------------------------------------------------------------------------ per file
def one_file(args):
    """Selection (v3 mc_select_file) plus universe weights of the selected events and truth sums."""
    path, unisims = args
    try:
        sel, truth = SEL.mc_select_file(path)
        root = uproot.open(path)
        S = root['phaseIITriggerTree']
    except Exception as e:                     # unreadable files are skipped, as in the C++ prep
        return dict(path=path, ok=False, why=f'{type(e).__name__}: {e}')
    with root:
        tm = SEL.truth_masks(S)
        w_ev, cv_raw = SEL.cv_weight(S)
        a0 = RW._jag(S, 'weight_All0_UBGenie', 100)
        kn = np.stack([RW._jag(S, f'weight_{k}', 100) for k in RW.FLUX_KNOBS])
        un = np.stack([RW._jag(S, f'weight_{s}_UBGenie', 2) for s in RW.UNISIM_SLOTS.values()]) if unisims else None
    cvr = np.where(np.isfinite(cv_raw), cv_raw, np.nan)
    flux = RW.safe(np.prod(kn, axis=0) * cvr[:, None])          # same product order as lib/reweight.py
    genie = RW.safe(a0)
    sig = np.where(tm['fv'] & tm['nu'])[0]
    out = dict(path=path, ok=True, truth=truth,
               T=dict(cv=float(w_ev[sig].sum()), x=genie[sig].sum(0), f=flux[sig].sum(0)))
    if unisims:
        out['T']['u'] = RW.safe(un[:, sig]).sum(1)            # (8, 2)
    if sel is not None:
        ce = sel['ce']
        sel['genie'], sel['flux'] = genie[ce], flux[ce]
        if unisims:
            sel['unisim'] = RW.safe(un[:, ce])                   # (8, n, 2)
        sel['fnum'] = np.full(len(ce), int(os.path.basename(path).split('_')[1].split('.')[0]))
        out['sel'] = sel
    return out


def run_files(files, unisims, workers):
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(one_file, [(p, unisims) for p in files], chunksize=8))
    bad = [r for r in res if not r['ok']]
    good = [r for r in res if r['ok']]
    sels = [r['sel'] for r in good if 'sel' in r]
    keys = [k for k in sels[0] if k != 'unisim']
    sel = {k: np.concatenate([s[k] for s in sels]) for k in keys}
    if unisims:
        sel['unisim'] = np.concatenate([s['unisim'] for s in sels], axis=1)
    T = dict(cv=sum(r['T']['cv'] for r in good), x=sum(r['T']['x'] for r in good),
             f=sum(r['T']['f'] for r in good))
    if unisims:
        T['u'] = sum(r['T']['u'] for r in good)
    truth = np.sum([r['truth'] for r in good], axis=0)
    return sel, T, truth, bad


# ------------------------------------------------------------------------------------------------ select
def cmd_select(a):
    files = [os.path.join(a.mc_dir, f) for f in read_list(a.list)] if a.mc_dir else read_list(a.list)
    sel, T, truth, bad = run_files(files, False, a.workers)
    # Per selected cluster (any bunch): (in window, in window & first 26 bunches, python category, OOFV nuNCQE)
    py_cmp = {(int(f), int(c)): (bool(w), bool(w and i), int(k), bool(o))
              for f, c, w, i, k, o in zip(sel['fnum'], sel['ce'], sel['in_win'], sel['in26'], sel['cat'],
                                          sel['oofv_nu'])}
    t = uproot.open(a.prep)['phaseIITriggerTree'].arrays(
        ['mc_file_number', 'trigger_entry', 'sel_ncqe_nobunch', 'sel_ncqe_inwin', 'sel_ncqe', 'ncqe_py_category',
         'category', 'ncqe_true_signal'], library='np')
    s = t['sel_ncqe_nobunch']
    cc_cmp = {(int(f), int(c)): (bool(w), bool(cv), int(k), int(cat) == 13)
              for f, c, w, cv, k, cat in zip(t['mc_file_number'][s], t['trigger_entry'][s], t['sel_ncqe_inwin'][s],
                                             t['sel_ncqe'][s], t['ncqe_py_category'][s], t['category'][s])}
    only_py, only_cc = set(py_cmp) - set(cc_cmp), set(cc_cmp) - set(py_cmp)
    diff = [k for k in set(py_cmp) & set(cc_cmp) if py_cmp[k] != cc_cmp[k]]
    n_sig_cc = int(t['ncqe_true_signal'].sum())
    print(f'files {len(files)} (unreadable {len(bad)})')
    print(f'selected clusters (all bunches): python {len(py_cmp)}  C++ {len(cc_cmp)}')
    print(f'  in window: python {sum(v[0] for v in py_cmp.values())}  C++ {sum(v[0] for v in cc_cmp.values())}')
    print(f'  CV (window & first 26): python {sum(v[1] for v in py_cmp.values())}  C++ {sum(v[1] for v in cc_cmp.values())}')
    cats = lambda d: {c: sum(1 for v in d.values() if v[1] and v[2] == c) for c in range(-1, 6)}
    print(f'  CV per category (py cat -1..5): python {cats(py_cmp)}  C++ {cats(cc_cmp)}')
    print(f'  CV OOFV nuNCQE: python {sum(v[1] and v[3] for v in py_cmp.values())}  '
          f'C++ {sum(v[1] and v[3] for v in cc_cmp.values())}')
    print(f'true FV nuNCQE triggers: python {int(truth[1])}  C++ {n_sig_cc}')
    print(f'only in python: {len(only_py)}  only in C++: {len(only_cc)}  flag differences: {len(diff)}')
    for k in sorted(only_py)[:10]:
        print('   only python', k, py_cmp[k])
    for k in sorted(only_cc)[:10]:
        print('   only C++   ', k, cc_cmp[k])
    for k in sorted(diff)[:10]:
        print('   differs    ', k, 'python', py_cmp[k], 'C++', cc_cmp[k])
    ok = not only_py and not only_cc and not diff and int(truth[1]) == n_sig_cc
    print('SELECTION MATCH' if ok else 'SELECTION MISMATCH')
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------ reference
def summarize(shifts):
    s = np.asarray(shifts) * 100
    d = dict(n=int(len(s)), rms=float(np.sqrt(np.mean(s ** 2))), mean=float(s.mean()), std=float(s.std()),
             shifts=[float(x) for x in s])
    if len(s) > 2:
        d.update(p16=float(np.percentile(s, 16)), p84=float(np.percentile(s, 84)))
    return d


def cmd_reference(a):
    files = [os.path.join(a.mc_dir, f) for f in read_list(a.list)] if a.mc_dir else read_list(a.list)
    sel, T, truth, bad = run_files(files, a.unisims, a.workers)
    N, D, Dsky = INPUTS['n_obs'], INPUTS['d_cit'], INPUTS['d_sky']
    fE = INPUTS['pot_data'] / a.mc_pot
    pot_1e20 = INPUTS['pot_data'] / 1e20

    Y0 = X.yields(sel, fE)
    f_cv = X.extract_f(Y0, N, D, Dsky)
    TE = truth[0] * fE
    sg = X.sigma_genie(TE, pot_1e20)
    sigma = f_cv * sg
    nb = N - D - Dsky - X.m_bkg(Y0)
    M = Y0['sig'] + Y0['OOFV_nu']

    def shifts(w, Tu, genie):
        Yu = X.yields(sel, fE, w=w)
        Yu['Ext'] = np.full_like(Yu['sig'], Y0['Ext'])          # External held at alpha x CV
        fu = X.extract_f(Yu, N, D, Dsky)
        return fu * ((Tu / T['cv']) if genie else 1.0) / f_cv - 1

    universes = {'weight_All_UBGenie': dict(summarize(shifts(sel['genie'], T['x'], True)), flux=False),
                 'weight_flux_all': dict(summarize(shifts(sel['flux'], T['f'], False)), flux=True)}
    if a.unisims:
        for j, slot in enumerate(RW.UNISIM_SLOTS.values()):
            universes[f'weight_{slot}_UBGenie'] = dict(summarize(shifts(sel['unisim'][j], T['u'][j], True)), flux=False)

    # MC statistics: delta method of 08_systematics.py
    cv = SEL.cv_mask(sel)
    w = np.where(np.isfinite(sel['w']), sel['w'], 0.0)
    cat, onu = sel['cat'], sel['oofv_nu']
    var = lambda mm, sc=1.0: (sc * fE) ** 2 * (w[cv & mm] ** 2).sum()
    vB = sum(var(mm) for mm in (cat == 1, cat == 2, cat == 3, (cat == 4) & ~onu)) + var(cat == 5, X.ALPHA_SCALING)
    vU = fE ** 2 * truth[5] - var(cat == 0)
    terms = dict(backgrounds=np.sqrt(vB) / nb, FV_signal=abs(1 / TE - 1 / M) * np.sqrt(var(cat == 0)),
                 OOFV_nuNCQE=np.sqrt(var(onu)) / M, truth_denominator=np.sqrt(vU) / TE)
    mc_tot = float(np.sqrt(sum(v ** 2 for v in terms.values())))

    raw = {c: int((cv & m).sum()) for c, m in dict(sig=cat == 0, OOFV_nu=onu, NCother=cat == 1, CC=cat == 2,
                                                  nubar=cat == 3, OOFV_other=(cat == 4) & ~onu, Ext=cat == 5).items()}
    out = dict(
        inputs=dict(files=len(files), unreadable=[b['path'] for b in bad], mc_pot=a.mc_pot, **INPUTS),
        cv=dict(yields={LABELS[k]: float(v) for k, v in Y0.items()}, raw_counts={LABELS[k]: v for k, v in raw.items()},
                M_bkg=float(X.m_bkg(Y0)), signal_like=float(M), numerator=float(nb), T=float(TE),
                sigma_genie=float(sg), f=float(f_cv), sigma=float(sigma)),
        universes=universes,
        mc_stat_terms_percent=dict({k: float(100 * v) for k, v in terms.items()}, total=100 * mc_tot),
        data_driven_percent=dict(stat=float(100 * np.sqrt(N) / nb), cit_stat=float(100 * INPUTS['d_cit_err'] / nb),
                                 skyshine=float(100 * INPUTS['d_sky_err'] / nb)),
        provenance=X.provenance(dict(script=os.path.abspath(__file__), v3=V3)))
    with open(a.out, 'w') as f:
        json.dump(out, f, indent=1)
    c = out['cv']
    print(f"files {len(files)} (unreadable {len(bad)}); raw CV counts {raw}")
    print(f"f = {c['f']:.6f}  sigma = {c['sigma']:.6f}  T = {c['T']:.4f}  M_bkg = {c['M_bkg']:.4f}")
    for k, v in universes.items():
        extra = f"p16 {v['p16']:.3f} p84 {v['p84']:.3f}" if 'p16' in v else f"shifts {v['shifts']}"
        print(f"  {k:40s} RMS {v['rms']:.4f}  mean {v['mean']:.4f}  {extra}")
    print(f"  MC stat {100 * mc_tot:.4f}%  ({', '.join(f'{k} {100 * v:.3f}' for k, v in terms.items())})")
    print(f"wrote {a.out}")
    return 0


# ------------------------------------------------------------------------------------------------ compare
def cmd_compare(a):
    cpp, ref = json.load(open(a.cpp)), json.load(open(a.ref))
    rows, worst = [], 0.0

    def add(name, x, y, tol):
        nonlocal worst
        rel = abs(x - y) / max(abs(y), 1e-300)
        worst = max(worst, rel) if tol else worst
        rows.append((name, x, y, rel, 'ok' if rel <= tol else 'DIFF'))

    def add_abs(name, x, y, tol):
        # absolute comparison, for quantities that can be zero (mean shifts of symmetric unisims)
        d = abs(x - y)
        rows.append((name, x, y, d, 'ok' if d <= tol else 'DIFF'))

    for k in ('f', 'sigma', 'T', 'M_bkg', 'signal_like', 'numerator', 'sigma_genie'):
        add(f'cv.{k}', cpp['cv'][k], ref['cv'][k], a.tol)
    for lab, v in ref['cv']['yields'].items():
        add(f'yield {lab}', cpp['cv']['yields'][lab], v, a.tol)
    for name, r in ref['universes'].items():
        if name not in cpp['universes']:
            rows.append((f'{name}', float('nan'), r['rms'], float('nan'), 'MISSING'))
            continue
        c = cpp['universes'][name]
        add(f'{name} RMS', c['rms'], r['rms'], a.tol)
        add_abs(f'{name} mean (abs diff, %)', c['mean'], r['mean'], a.tol_abs)
        shifts_c, shifts_r = np.array(c['shifts']), np.array(r['shifts'])
        rows.append((f'{name} max|shift diff| (%)', float(np.max(np.abs(shifts_c - shifts_r))), 0.0, float('nan'), ''))
    add('MC stat total', cpp['mc_stat_terms_percent']['total'], ref['mc_stat_terms_percent']['total'], a.tol)
    for k in ('stat', 'cit_stat', 'skyshine'):
        add(f'{k}', cpp['data_driven_percent'][k], ref['data_driven_percent'][k], a.tol)
    print(f"{'quantity':52s} {'C++':>16s} {'python':>16s} {'diff':>10s}  (relative unless marked abs)")
    for n, x, y, r, s in rows:
        print(f'{n:52s} {x:16.8g} {y:16.8g} {r:10.2e} {s}')
    ok = all(s != 'DIFF' and s != 'MISSING' for *_, s in rows)
    print('NUMBERS MATCH' if ok else 'NUMBERS DIFFER')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('select', 'reference'):
        p = sub.add_parser(name)
        p.add_argument('--mc-dir', default='', help='directory prefixed to the names in --list')
        p.add_argument('--list', required=True, help='file names (one per line), as given to annie_ncqe_prep')
        p.add_argument('--workers', type=int, default=int(os.environ.get('NCQE_N_WORKERS', '8')))
    sub.choices['select'].add_argument('--prep', required=True, help='merged annie_ncqe_prep output')
    r = sub.choices['reference']
    r.add_argument('--mc-pot', type=float, required=True, help='simulated POT of the listed files')
    r.add_argument('--unisims', action='store_true')
    r.add_argument('--out', required=True)
    c = sub.add_parser('compare')
    c.add_argument('cpp')
    c.add_argument('ref')
    c.add_argument('--tol', type=float, default=1e-9, help='relative tolerance')
    c.add_argument('--tol-abs', type=float, default=1e-9, help='absolute tolerance on mean shifts (percent)')
    a = ap.parse_args()
    return dict(select=cmd_select, reference=cmd_reference, compare=cmd_compare)[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
