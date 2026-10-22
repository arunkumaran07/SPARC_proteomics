#' Supported species
#'
#' Species determine the annotation package used to map gene symbols to Entrez
#' IDs and the organism name passed to Reactome. Choose `"other"` for species
#' without Reactome support; enrichment is then skipped and everything else
#' runs normally.
#'
#' @return A data frame with one row per species: the config key, the
#'   Bioconductor `OrgDb` package, the Reactome organism name used by
#'   \pkg{ReactomePA}, and a display name.
#' @export
#' @examples
#' list_species()
list_species <- function() {
  data.frame(
    species  = c("human", "mouse", "rat", "zebrafish", "fly", "worm", "other"),
    orgdb    = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
                 "org.Dr.eg.db", "org.Dm.eg.db", "org.Ce.eg.db", NA),
    reactome = c("human", "mouse", "rat", "zebrafish", "fly", "celegans", NA),
    label    = c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                 "Danio rerio", "Drosophila melanogaster",
                 "Caenorhabditis elegans", "Other / not annotated"),
    stringsAsFactors = FALSE
  )
}

#' Look up one species
#'
#' @param species A key from [list_species()]. Common aliases such as
#'   `"Homo sapiens"`, `"hs"`, `"Mus musculus"` or `"mm"` are accepted.
#' @return A one-row data frame from [list_species()].
#' @export
#' @examples
#' species_info("Mus musculus")
species_info <- function(species) {
  aliases <- c("homo sapiens" = "human", hs = "human", hsa = "human",
               "mus musculus" = "mouse", mm = "mouse", mmu = "mouse",
               "rattus norvegicus" = "rat", rn = "rat", rno = "rat",
               "danio rerio" = "zebrafish", dr = "zebrafish",
               "drosophila melanogaster" = "fly", dm = "fly",
               "caenorhabditis elegans" = "worm", celegans = "worm", ce = "worm")
  key <- tolower(trimws(species))
  if (key %in% names(aliases)) key <- aliases[[key]]
  tab <- list_species()
  if (!key %in% tab$species)
    stop("Unknown species '", species, "'. Choose one of: ",
         paste(tab$species, collapse = ", "), call. = FALSE)
  tab[tab$species == key, , drop = FALSE]
}
