#' certusR: likelihood-ratio parentage inference
#'
#' A dependency-free implementation of likelihood-ratio parentage inference for
#' codominant markers, reconstructed from the primary literature.
#'
#' @section The name:
#' Latin \emph{certus}, "determined, resolved, certain": the past participle of
#' \emph{cernere}, to sift apart, to distinguish, to decide. The same root gives
#' English \emph{discern}, \emph{criterion} and \emph{crisis}. It is also one
#' letter from \emph{cervus}, the red deer of the Isle of Rum whose pedigree
#' motivated Marshall et al. (1998) and gave CERVUS its name. The distinction is
#' the point: CERVUS returns a score, and what this package adds is how far that
#' score can be trusted for the particular trio in front of you.
#'
#' @section Workflow:
#' \enumerate{
#'   \item [genotypes()] to read a wide allele table.
#'   \item [allele_freqs()] for reference frequencies.
#'   \item [estimate_error_rates()] for a per-locus genotyping error rate,
#'         if known parent-offspring pairs are available.
#'   \item [parentage_lod()] for LOD scores, and [delta_stat()] for Delta.
#'   \item [assign_parentage()] to attach confidence, using either the
#'         trio-specific criterion of [lod_null()] and [lod_critical()] or the
#'         population-level criterion of [sim_delta()] and [delta_critical()].
#' }
#'
#' @section Provenance:
#' The likelihood equations were reconstructed from the published equations and
#' prose of the three papers cited below. No existing implementation's source
#' code was consulted. The error expansion is re-derived from the stated
#' random-genotype-replacement model rather than transcribed, because the
#' coefficients printed in the sources do not form a probability decomposition;
#' a corrigendum to Kalinowski et al. (2007) records typesetting errors in that
#' paper's appendix, and the coefficients of Amiri Roudbar et al. (2025) show the
#' same damage. The full argument, with the term-by-term derivation, is in
#' `system.file("notes", "derivation.md", package = "certusR")`.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x}
#'
#' Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
#' (2025) Pairwise paternity assignment with forward-backward simulations:
#' refining CERVUS using trio-based likelihood and locus-specific error rates.
#' \emph{Ecology and Evolution}, \strong{15}(10), e72230.
#' \doi{10.1002/ece3.72230}
#'
#' @importFrom stats median quantile runif
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"
