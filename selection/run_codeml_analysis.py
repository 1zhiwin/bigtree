#!/usr/bin/env python3
"""
Run CodeML (PAML) positive selection analysis on DAH7PS sequences.
"""

import os
import subprocess
import re
from pathlib import Path
import json

def create_codeml_control_file(alignment_file, tree_file, output_prefix, model_params):
    """
    Create CodeML control file for a specific model.

    model_params: dict with 'NSsites', 'model', 'fix_omega', 'omega'
    """

    control_content = f"""      seqfile = {alignment_file}
     treefile = {tree_file}
      outfile = {output_prefix}.out

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 2    * 2: amino acids
    CodonFreq = 2    * 2: F3X4

        model = {model_params.get('model', 0)}
      NSsites = {model_params.get('NSsites', 0)}

        icode = 0
        Mgene = 0

    fix_kappa = 0
        kappa = 2
    fix_omega = {model_params.get('fix_omega', 0)}
        omega = {model_params.get('omega', 0.4)}

    fix_alpha = 1
        alpha = 0.
       Malpha = 0
        ncatG = 8

        getSE = 0
 RateAncestor = 0

   Small_Diff = .5e-6
    cleandata = 1
       method = 0
"""

    control_file = f"{output_prefix}.ctl"
    with open(control_file, 'w') as f:
        f.write(control_content)

    return control_file

def run_codeml(control_file, work_dir="."):
    """Run CodeML with the specified control file."""

    original_dir = os.getcwd()
    os.chdir(work_dir)

    try:
        # Run codeml
        result = subprocess.run(['codeml', os.path.basename(control_file)],
                              capture_output=True,
                              text=True,
                              timeout=600)  # 10 minute timeout

        return_code = result.returncode
        stdout = result.stdout
        stderr = result.stderr

    except subprocess.TimeoutExpired:
        print(f"  WARNING: CodeML timed out after 10 minutes")
        return_code = -1
        stdout = ""
        stderr = "Timeout"
    finally:
        os.chdir(original_dir)

    return return_code, stdout, stderr

def parse_codeml_output(output_file):
    """Parse CodeML output file to extract key results."""

    if not os.path.exists(output_file):
        print(f"  ERROR: Output file {output_file} not found")
        return None

    with open(output_file, 'r') as f:
        content = f.read()

    results = {}

    # Extract lnL (log-likelihood)
    lnl_match = re.search(r'lnL\(ntime:\s+\d+\s+np:\s+\d+\):\s+([-\d.]+)', content)
    if lnl_match:
        results['lnL'] = float(lnl_match.group(1))

    # Extract omega (dN/dS)
    omega_match = re.search(r'omega \(dN/dS\) =\s+([\d.]+)', content)
    if omega_match:
        results['omega'] = float(omega_match.group(1))

    # Extract number of parameters
    np_match = re.search(r'lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)', content)
    if np_match:
        results['np'] = int(np_match.group(1))

    # Extract site class parameters (for M2a, M8)
    # Look for proportion and omega values
    prop_match = re.findall(r'proportion\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', content)
    if prop_match:
        results['proportions'] = prop_match

    # Extract positively selected sites (BEB analysis)
    positively_selected = []
    beb_section = re.search(r'Bayes Empirical Bayes \(BEB\) analysis(.*?)(\n\n|\nTime used)', content, re.DOTALL)
    if beb_section:
        beb_text = beb_section.group(1)
        # Parse sites with Pr(w>1) > 0.95
        site_lines = re.findall(r'(\d+)\s+([A-Z])\s+([\d.]+)\s+\+\-\s+([\d.]+)\s+([\*\s]+)', beb_text)
        for site, aa, omega, se, sig in site_lines:
            if '*' in sig or '**' in sig:  # * or ** indicates significant
                positively_selected.append({
                    'site': int(site),
                    'amino_acid': aa,
                    'omega': float(omega),
                    'se': float(se),
                    'significance': sig.count('*')
                })

    results['positively_selected_sites'] = positively_selected

    return results

def run_site_models(alignment_file, tree_file, output_dir, type_name):
    """Run site models M0, M1a, M2a, M7, M8."""

    print(f"\n{'='*60}")
    print(f"Running Site Models: {type_name}")
    print(f"{'='*60}\n")

    models = {
        'M0': {'NSsites': 0, 'model': 0},  # One-ratio
        'M1a': {'NSsites': 1, 'model': 0},  # Nearly neutral
        'M2a': {'NSsites': 2, 'model': 0},  # Positive selection
        'M7': {'NSsites': 7, 'model': 0},  # Beta distribution
        'M8': {'NSsites': 8, 'model': 0},  # Beta + omega
    }

    results = {}

    for model_name, params in models.items():
        print(f"Running {model_name}...")

        output_prefix = os.path.join(output_dir, model_name)
        control_file = create_codeml_control_file(
            alignment_file, tree_file, output_prefix, params
        )

        # Run CodeML
        return_code, stdout, stderr = run_codeml(control_file, output_dir)

        if return_code == 0:
            # Parse output
            output_file = f"{output_prefix}.out"
            model_results = parse_codeml_output(output_file)

            if model_results:
                results[model_name] = model_results
                print(f"  {model_name} completed successfully")
                print(f"    lnL = {model_results.get('lnL', 'N/A')}")
                if 'omega' in model_results:
                    print(f"    omega = {model_results['omega']:.4f}")
                if model_results.get('positively_selected_sites'):
                    print(f"    Positively selected sites: {len(model_results['positively_selected_sites'])}")
            else:
                print(f"  {model_name} completed but could not parse output")
        else:
            print(f"  {model_name} failed with return code {return_code}")
            if stderr:
                print(f"  Error: {stderr[:200]}")

    return results

def perform_likelihood_ratio_tests(results):
    """Perform likelihood ratio tests between nested models."""

    print(f"\n{'='*60}")
    print(f"Likelihood Ratio Tests")
    print(f"{'='*60}\n")

    lrt_results = {}

    # Test M1a vs M2a
    if 'M1a' in results and 'M2a' in results:
        lnL1 = results['M1a']['lnL']
        lnL2 = results['M2a']['lnL']
        df = results['M2a']['np'] - results['M1a']['np']
        lrt = 2 * (lnL2 - lnL1)

        # Chi-square critical value at p=0.05, df=2 is ~5.99
        # Chi-square critical value at p=0.01, df=2 is ~9.21
        significant = "Yes" if lrt > 5.99 else "No"

        print(f"M1a vs M2a (Positive Selection Test):")
        print(f"  M1a lnL: {lnL1:.2f}")
        print(f"  M2a lnL: {lnL2:.2f}")
        print(f"  LRT statistic: {lrt:.4f} (df={df})")
        print(f"  Significant (p<0.05): {significant}")

        lrt_results['M1a_vs_M2a'] = {
            'LRT': lrt,
            'df': df,
            'significant_p05': lrt > 5.99,
            'significant_p01': lrt > 9.21
        }

    # Test M7 vs M8
    if 'M7' in results and 'M8' in results:
        lnL7 = results['M7']['lnL']
        lnL8 = results['M8']['lnL']
        df = results['M8']['np'] - results['M7']['np']
        lrt = 2 * (lnL8 - lnL7)

        significant = "Yes" if lrt > 5.99 else "No"

        print(f"\nM7 vs M8 (Positive Selection Test):")
        print(f"  M7 lnL: {lnL7:.2f}")
        print(f"  M8 lnL: {lnL8:.2f}")
        print(f"  LRT statistic: {lrt:.4f} (df={df})")
        print(f"  Significant (p<0.05): {significant}")

        lrt_results['M7_vs_M8'] = {
            'LRT': lrt,
            'df': df,
            'significant_p05': lrt > 5.99,
            'significant_p01': lrt > 9.21
        }

    return lrt_results

if __name__ == "__main__":

    print("\n" + "="*60)
    print("CODEML POSITIVE SELECTION ANALYSIS")
    print("="*60)

    # Type I analysis
    type_i_alignment = "../msa/type_i/type_i_alignment_trimmed.faa"
    type_i_tree = "../trees/type_i/type_i_tree.treefile"
    type_i_output = "type_i"

    type_i_results = run_site_models(
        type_i_alignment,
        type_i_tree,
        type_i_output,
        "Type I"
    )

    # Save Type I results
    with open("type_i/codeml_results.json", 'w') as f:
        # Convert to JSON-serializable format
        json_results = {}
        for model, data in type_i_results.items():
            json_results[model] = {
                k: v for k, v in data.items()
                if k != 'proportions'  # Skip complex nested structures
            }
        json.dump(json_results, f, indent=2)

    # LRT tests for Type I
    type_i_lrt = perform_likelihood_ratio_tests(type_i_results)

    with open("type_i/lrt_results.json", 'w') as f:
        json.dump(type_i_lrt, f, indent=2)

    # Type II analysis
    print("\n")
    type_ii_alignment = "../msa/type_ii/type_ii_alignment_trimmed.faa"
    type_ii_tree = "../trees/type_ii/type_ii_tree.treefile"
    type_ii_output = "type_ii"

    type_ii_results = run_site_models(
        type_ii_alignment,
        type_ii_tree,
        type_ii_output,
        "Type II"
    )

    # Save Type II results
    with open("type_ii/codeml_results.json", 'w') as f:
        json_results = {}
        for model, data in type_ii_results.items():
            json_results[model] = {
                k: v for k, v in data.items()
                if k != 'proportions'
            }
        json.dump(json_results, f, indent=2)

    # LRT tests for Type II
    type_ii_lrt = perform_likelihood_ratio_tests(type_ii_results)

    with open("type_ii/lrt_results.json", 'w') as f:
        json.dump(type_ii_lrt, f, indent=2)

    print("\n" + "="*60)
    print("CODEML ANALYSIS COMPLETE")
    print("="*60)
    print("\nResults saved to:")
    print("  - selection/type_i/codeml_results.json")
    print("  - selection/type_i/lrt_results.json")
    print("  - selection/type_ii/codeml_results.json")
    print("  - selection/type_ii/lrt_results.json")
