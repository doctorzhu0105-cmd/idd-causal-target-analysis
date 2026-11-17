"""Portable entry point for the recovered, bounded analysis components.

This is an engineering dispatcher. Numerical routines live in the recovered
modules. Full catalogue and historical source replays are not implied.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
sys.path.insert(0, str(HERE))

SUPPORTED = ('givw', 'attenuation', 'k1-propagation', 'k-transition', 'r39-primary',
             'qwen-ab', 'qwen-rho', 'qwen-r37', 'qwen-psd',
             'figure-1', 'figure-2', 'figure-s1', 'figure-s2', 'figure-s3')

def load_file(directory: Path, name: str) -> Path:
    if not isinstance(name, str) or not name:
        raise ValueError('A relative input file name is required in the config')
    path = (directory / name).resolve()
    if not path.is_relative_to(directory) or not path.is_file():
        raise ValueError('Required input file is missing or outside input-dir: ' + name)
    return path

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', choices=SUPPORTED, required=True)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    try:
        inputs, out, config_file = args.input_dir.resolve(), args.output_dir.resolve(), args.config.resolve()
        if not inputs.is_dir() or not config_file.is_file():
            raise ValueError('input-dir or config does not exist')
        if out.exists():
            raise ValueError('output-dir must not exist')
        if out.is_relative_to(inputs) or inputs.is_relative_to(out):
            raise ValueError('output-dir must be disjoint from input-dir')
        if out.is_relative_to(PACKAGE) or PACKAGE.is_relative_to(out):
            raise ValueError('output-dir must be outside the release package')
        if config_file.is_relative_to(out):
            raise ValueError('config cannot be inside the output directory')
        config = json.loads(config_file.read_text(encoding='utf-8-sig'))
        if not isinstance(config, dict):
            raise ValueError('config must be a JSON object')

        if args.analysis in ('givw', 'attenuation'):
            import numpy as np
            from project_modules.correlated_mr import ivw_correlated
            from run_scientific_analysis import attenuation
            source = load_file(inputs, config.get('input'))
            data = json.loads(source.read_text(encoding='utf-8-sig'))
            if args.analysis == 'attenuation':
                result = attenuation(data, config.get('candidate'), float(config['lambda']))
            else:
                r = ivw_correlated(np.asarray(data['beta_x']), np.asarray(data['beta_y']),
                                   np.asarray(data['se_y']), np.asarray(data['signed_ld']))
                result = {k: (None if isinstance(v, float) and not math.isfinite(v) else v)
                          for k, v in r.as_dict().items()}
            out.mkdir(parents=True, exist_ok=False)
            (out/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        elif args.analysis == 'k1-propagation':
            import pandas as pd
            from run_scientific_analysis import propagate
            fdr = load_file(inputs, config.get('fdr'))
            official = load_file(inputs, config.get('official'))
            calculated = propagate(pd.read_csv(fdr, sep='\t', comment='#'),
                                   pd.read_csv(official, sep='\t', comment='#'))
            out.mkdir(parents=True, exist_ok=False)
            calculated.to_csv(out/'k1_second_order_multiplicity.tsv', sep='\t', index=False)
            result = {'rows':len(calculated), 'k1_rows':int((calculated.k == 1).sum()),
                      'bh_changes':int(calculated.bh_membership_changed.sum()),
                      'by_changes':int(calculated.by_membership_changed.sum()),
                      'scope':'global adjustment of supplied estimable rows; not upstream MR generation'}
            (out/'result.json').write_text(json.dumps(result, indent=2)+'\n',encoding='utf-8')
        elif args.analysis == 'r39-primary':
            from r39_primary import compute
            source = load_file(inputs, config.get('input'))
            result = compute(source)
            out.mkdir(parents=True, exist_ok=False)
            (out/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        elif args.analysis in ('qwen-ab', 'qwen-rho', 'qwen-r37', 'qwen-psd'):
            # Independent saved-result audits, not original analysis generators.
            if args.analysis == 'qwen-ab':
                from verification.qwen_ab import run
                first = load_file(inputs, config.get('paired'))
                second = load_file(inputs, config.get('reference'))
                out.mkdir(parents=True, exist_ok=False)
                run(first, second, out)
            else:
                if args.analysis == 'qwen-rho':
                    from verification.qwen_rho import run
                elif args.analysis == 'qwen-r37':
                    from verification.qwen_r37 import run
                else:
                    from verification.qwen_psd import run
                first = load_file(inputs, config.get('input'))
                out.mkdir(parents=True, exist_ok=False)
                run(first, out)
            result = {'analysis': args.analysis, 'scope': 'independent numerical verification of supplied saved inputs'}
        elif args.analysis == 'k-transition':
            from run_release_analysis import k_transition
            source = load_file(inputs, config.get('input'))
            out.mkdir(parents=True, exist_ok=False)
            result = k_transition(source, out)
        else:
            spec = importlib.util.spec_from_file_location('historical_v5_figures', HERE/'figures'/'v5_figures.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            keys = {'figure-1':['effective_k','exposure','precision'], 'figure-2':['figure2'],
                    'figure-s1':[], 'figure-s2':['figure_s2'], 'figure-s3':['leverage']}[args.analysis]
            files = {k:load_file(inputs, config.get(k)) for k in keys}
            out.mkdir(parents=True, exist_ok=False)
            module.generate(args.analysis, files, out)
            result = {'figure':args.analysis, 'source':'adapted v5 figure code',
                      'notice':'font and rendering environment can change output bytes'}
            (out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(result, ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        print('INPUT_OR_EXECUTION_ERROR: '+str(exc), file=sys.stderr)
        return 2

if __name__ == '__main__':
    raise SystemExit(main())
