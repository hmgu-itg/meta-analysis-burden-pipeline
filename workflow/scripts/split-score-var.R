library(data.table)

options(warn = 2)
# This is so that a warning message is instead raised as an
# error and terminates the script when encountered.

score.files = snakemake@input[['scores']]
var.files = snakemake@input[['vars']]
group.file = snakemake@input[['groupfile']]

chromosome = snakemake@wildcards[['chrom']]


n.files = length(score.files)


group.info = fread(group.file, header = FALSE, colClasses = list(
        character = c(1, 2, 4, 5),
        numeric = c(3, 6)
    ))
setnames(group.info, c("group", "chr", "pos", "ref", "alt", "weight"))

group.info = group.info[chr==chromosome]
chrom.groups = unique(group.info$group)

group.info[,variant.id:=paste(group, chr, pos, ref, alt, sep = ":")]
if (any(duplicated(group.info$variant.id))) {
    warning("Duplicated variant ID detected in group.file. Excluding all duplicates, but this should be checked!!")
    group.info <- group.info[!duplicated(variant.id), ]
}

# Below is needed for memory efficient subsetting of the group.info later using variant.ids
setkey(group.info, variant.id)


cohort.scores <- NULL
cohort.covms = list()
for(i in 1:n.files) {
    score.file = score.files[i]
    score <- fread(
        score.file, 
        header = TRUE, 
        colClasses = list(
            character = c(1, 4, 5), 
            integer = c(2, 3, 6), 
            numeric = 7:11
        )
    )
    score <- score[,.(group, chr, pos, ref, alt, N, missrate, altfreq, SCORE)]
    score[,variant.id:=paste(group, chr, pos, ref, alt, sep = ":")]
    cohort.scores <- rbind(cohort.scores, score)

    if (!any(score$variant.id %in% group.info$variant.id)) {
        missing.variant.count = table(score$variant.id %in% group.info$variant.id)[["FALSE"]]
        cat(missing.variant.count,
                "variants in the",
                score.file,
                "are missing in the group file.")
        stop("Error: meta files possibly not generated using this group.file!")
    }

    var.file = paste0(var.files[i])
    con = file(var.file, "rb")
    for (selected.group in unique(score$group)) {
        # `unique` does not inherently sort the unique list. The ordering is preserved,
        # so it should be okay to just simply loop through instead of indexing
        group.n.p = nrow(score[group==selected.group])
        group.V = matrix(0, group.n.p, group.n.p)
        group.V[lower.tri(group.V, diag = TRUE)] <- readBin(con, what = "numeric", n = (1+group.n.p)*group.n.p/2, size = 4)
        # `options(warn = 2)` is set for the above line
        # so that it terminates when the var.file runs out of values to read
        # while there is still values to be filled in the grou.V matrix (Read the error.md for more details).
        group.V[upper.tri(group.V)] <- t(group.V)[upper.tri(group.V)]
        if (all(group.V == 0)) {
            stop(paste("Failed to recreate covariance matrix for", group, "in", var.file))
        }
        cohort.covms[[selected.group]] = group.V

    }
    remaining.check = readBin(con, what = "numeric", n = 1, size = 4)
    if (length(remaining.check)>0) {
        cat("Values still remaining in", var.file, "after recreating all covariance matrix.")
        stop(paste("Mismatching .score and .var file used:", score.file, "and", var.file))
    }
    close(con)
}

cohort.scores = cohort.scores[group %in% chrom.groups]
cohort.scores[,`:=`(
    Nmiss = (N*missrate/(1-missrate)),
    AC = (2*N*altfreq)
    )]

fwrite(cohort.scores, file = snakemake@output[['combined_scores']], quote = FALSE)
save(cohort.covms, file = snakemake@output[['combined_covms']])