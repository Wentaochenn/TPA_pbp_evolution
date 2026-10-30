from collections import defaultdict
from pathlib import Path
import math
from prepare import read, table
FOCAL = ['mrcA:p.Ala506Thr', 'mrcA:p.Ala506Val', 'mrcA:p.Gly708Ser', 'mrcA:p.Met625Val', 'pbp1:p.Pro564Leu', 'tp47:p.Ser394Arg']

def number(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return math.nan

def jaccard(a, b):
    x = set(filter(None, a.split(';')))
    y = set(filter(None, b.split(';')))
    return len(x & y)/len(x | y) if x | y else 0

def compare_events(work):
    work = Path(work)
    data = {model: read(work/'clock_sensitivity'/model/'pbp_dated_independent_origins.tsv')
            for model in ('strict','relaxed')}
    results = []
    for locus in FOCAL:
        a = [r for r in data['strict'] if r['locus'] == locus]
        b = [r for r in data['relaxed'] if r['locus'] == locus]
        used = set()
        for i, x in enumerate(a):
            scores = [jaccard(x['carrier_tip_set'],y['carrier_tip_set']) for y in b]
            best = max(scores,default=0)
            candidates = [j for j,s in enumerate(scores) if abs(s-best)<1e-12 and s>=.5]
            match = None
            status = 'no_matching_relaxed_carrier_set'
            if len(candidates)==1:
                j = candidates[0]
                reverse = [jaccard(y['carrier_tip_set'],b[j]['carrier_tip_set']) for y in a]
                reverse_best = max(reverse)
                if sum(abs(s-reverse_best)<1e-12 for s in reverse)==1 and reverse[i]==reverse_best and j not in used:
                    match=b[j];used.add(j);status='reciprocal_unique_match'
                else:
                    status='ambiguous_carrier_set_match'
            elif len(candidates)>1:
                status='ambiguous_carrier_set_match'
            row = {'locus':locus,'strict_event':x['event_uid'],
                   'relaxed_event':match['event_uid'] if match else '',
                   'carrier_set_jaccard':best,'matching_status':status}
            for model,r in [('strict',x),('relaxed',match)]:
                for key in ('posterior_simmap_median','posterior_simmap_hpd_low','posterior_simmap_hpd_high',
                            'posterior_simmap_p_after_1943','posterior_event_matching_support','posterior_n_assigned'):
                    row[model+'_'+key]=r.get(key,'NA') if r else 'NA'
            d1=number(row['strict_posterior_simmap_median']);d2=number(row['relaxed_posterior_simmap_median'])
            row['relaxed_minus_strict_years']=d2-d1
            results.append(row)
        for j,y in enumerate(b):
            if j not in used:
                results.append({'locus':locus,'strict_event':'','relaxed_event':y['event_uid'],
                                'matching_status':'relaxed_event_not_uniquely_matched_to_strict'})
    keys=list(dict.fromkeys(k for r in results for k in r))
    table(work/'clock_sensitivity/origin_clock_comparison.tsv',results,keys)

    probabilities=[]
    for model in ('strict','relaxed'):
        folder=work/'clock_sensitivity'/model
        maps=read(folder/'pbp_mapping_execution_audit.tsv')
        events=read(folder/'pbp_posterior_simmap_gain_events.tsv')
        by=defaultdict(list)
        for e in events:
            by[e['locus'],e['posterior_tree_index']].append(number(e['event_date']))
        for locus in FOCAL:
            successful=[m for m in maps if m['locus']==locus and m['status']=='complete']
            with_gain=[]
            for m in successful:
                dates=by[locus,m['posterior_tree_index']]
                if dates:
                    with_gain.append(min(dates)>1943)
            n=len(successful);g=len(with_gain)
            probabilities.append({'model':model,'locus':locus,'successful_tree_maps':n,
                'maps_with_internal_gain':g,'maps_without_internal_gain':n-g,
                'p_at_least_one_internal_gain':g/n if n else 'NA',
                'p_all_internal_gains_after_1943_given_any_gain':sum(with_gain)/g if g else 'NA',
                'p_any_gain_and_all_internal_gains_after_1943':sum(with_gain)/n if n else 'NA',
                'conditioning':'ER fitted on CA MCC; root prior 0=.999,1=.001; internal gains only'})
    table(work/'clock_sensitivity/locus_joint_probability.tsv',probabilities)
