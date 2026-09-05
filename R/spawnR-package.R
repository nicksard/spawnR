#' spawnR: likelihood-ratio parentage inference
#'
#' A dependency-free implementation of likelihood-ratio parentage inference for
#' codominant markers, reconstructed from the primary literature.
#'
#' @section The name:
#' To \emph{spawn}, of an aquatic animal, is to shed eggs or milt into open
#' water; from Latin \emph{expandere}, to spread out. Spawning is reproduction
#' that cannot be watched. Gametes are broadcast, several males may fertilise one
#' female's clutch at once, and no pair is visible from the bank. That
#' unobservability is the reason genetic parentage assignment exists in these
#' systems at all, and it is what this package is for.
#'
#' It reads the other way too: to \code{spawn} is what a program does when it
#' starts a child process. Both readings are true, which is the point.
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
