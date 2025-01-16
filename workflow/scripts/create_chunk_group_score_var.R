library(data.table)

output.score.file = snakemake@output[['score']]
output.var.file = snakemake@output[['var']]

cohort_ = snakemake@wildcards[['cohort']]
chrom_ = snakemake@wildcards[['chrom']]


subset.groups = fread(snakemake@input[['groups']], header=F)$V1

# Load, subset, and output score file
cohort.scores = fread(snakemake@input[['combined_scores']])

score = cohort.scores[
        chr == chrom_
        & group %in% subset.groups
    ]
rm(cohort.scores)


# Load, subset, and output var file
load(snakemake@input[['combined_covms']])
groups.to.include = names(cohort.covms) %in% subset.groups
chunk.covms = cohort.covms[groups.to.include]

cohort.data = list(
    score = score,
    covms = chunk.covms
)
save(cohort.data, file = snakemake@output[['data']])

load(snakemake@output[['data']]) # Try loading straight after saving to ensure it saved correctly