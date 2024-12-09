include: 'prepare-groups.smk'


def split_scorevar_input(w):
    input_prefix_series = run_list.loc[
                                (run_list['phenotype']==w.phenotype)
                                & (run_list['chrom'].astype(str)==str(w.chrom))
                                & (run_list['cohort']==w.cohort), 'input-prefix']
    if input_prefix_series.shape[0]!=1:
        raise ValueError(f"Unexpected number of input-prefix for {w.phenotype}, chr{w.chrom}, {w.cohort}")
    input_prefix = input_prefix_series.iloc[0]

    _directory = Path(input_prefix).parent
    file_prefix = Path(input_prefix).name
    scores = sorted(list(_directory.rglob(f"{file_prefix}.score.*")))
    _vars = sorted(list(_directory.rglob(f"{file_prefix}.var.*")))
    if len(scores) != len(_vars):
        raise ValueError(f"Differing number of score and var files for {w.phenotype}, chr{w.chrom}, {w.cohort}")
    return scores, _vars


rule split_scorevar:
    input:
        scores=lambda w: split_scorevar_input(w)[0],
        vars=lambda w: split_scorevar_input(w)[1],
        groupfile=config['group.file']
    output:
        combined_scores=temp('output/{phenotype}/per-cohort/{cohort}/chr{chrom}/all.scores.txt'),
        combined_covms=temp('output/{phenotype}/per-cohort/{cohort}/chr{chrom}/all.covms.RData')
    container: config['container']
    script:
        '../scripts/split-score-var.R'

rule create_chunk_score_var:
    input:
        groups=rules.create_chunk_groups.output[0],
        combined_scores=rules.split_scorevar.output.combined_scores,
        combined_covms=rules.split_scorevar.output.combined_covms
    output:
        data=temp('output/{phenotype}/per-cohort/{cohort}/chr{chrom}/chunk_{num}.RData')
    container: config['container']
    script:
        '../scripts/create_chunk_group_score_var.R'


def run_input(w):
    cohorts = run_list.loc[(run_list['phenotype']==w.phenotype)
                 & (run_list['chrom'].astype(str)==str(w.chrom)), 'cohort'].drop_duplicates()
    return expand(rules.create_chunk_score_var.output.data, cohort=cohorts, allow_missing=True)

rule run:
    input:
        data=run_input,
        groups=rules.create_chunk_groups.output[0],
        groupfile=config['group.file']
    params:
        SMMAT_config=config['SMMAT.config']
    threads: workflow.cores
    resources:
        mem_mb=config['resources']['run']['mem-mb']
    log:
        'output/{phenotype}/chr{chrom}/chunk_{num}.log'
    output:
        result=temp('output/{phenotype}/chr{chrom}/chunk_{num}.result')
    container: config['container']
    script:
        '../scripts/analyse.R'


rule debug_run_single_group:
    input:
        data=run_input,
        groups=rules.create_chunk_groups.output[0],
        groupfile=config['group.file']
    params:
        SMMAT_config=config['SMMAT.config']
    threads: workflow.cores
    resources:
        mem_mb=config['resources']['run']['mem-mb']
    output:
        result=temp('output/{phenotype}/debug/chr{chrom}/debug.chunk_{num}.{group}.result')
    container: config['container']
    script:
        '../scripts/analyse_single_group.R'


def collate_chrom_input(w):
    chunk_count = pd.read_csv(checkpoints.create_all_chunk_groups.get(chrom=w.chrom).output[0],
                                header = None)[0][0]
    chunk_nums = range(0, chunk_count)
    return expand(rules.run.output.result, num=chunk_nums, allow_missing=True)

rule collate_chrom:
    input:
        collate_chrom_input
    output:
        'output/{phenotype}/chr{chrom}.results.txt'
    run:
        combined = pd.concat([pd.read_csv(f, header = 0) for f in input])
        unique_groups = group_info.loc[group_info["chrom"].astype(str)==str(wildcards.chrom), 'group'].drop_duplicates()
        combined['group'] = pd.Categorical(combined['group'], categories=unique_groups, ordered=True)
        combined.sort_values('group', inplace=True)
        combined.to_csv(output[0], header = True, index = False)


def collate_phenotype_input(w):
    chroms = run_list.loc[run_list['phenotype']==w.phenotype, 'chrom'].drop_duplicates()
    return expand(rules.collate_chrom.output[0], chrom=chroms, allow_missing=True)

rule collate_phenotype:
    input:
        collate_phenotype_input
    output:
        'output/{phenotype}.csv'
    run:
        combined = pd.concat([pd.read_csv(f, header = 0) for f in input])
        unique_groups = run_list.loc[run_list["chrom"].astype(str)==str(wildcards.chrom), 'group'].drop_duplicates()
        combined['group'] = pd.Categorical(combined['group'], categories=unique_groups, ordered=True)
        combined.sort_values('group', inplace=True)
        combined.to_csv(output[0], header = True, index = False)


    