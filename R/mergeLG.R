##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#'
#' FS method: Merge Linkage Groups
#' 
#' Method for merging multiple linkage groups contained
#' in an FS object.
#' 
#' For a linkage analysis in GUSMap, there may be need to merge linkage groups. The indices of the 
#' linkage groups corresponds to the number
#' given in the output for the FS object. There needs to be at least two linkage groups 
#' specified for merging and more than two linkage groups can be merged at once.
#' 
#' When \code{where = "LG-pts"}, MI and PI linkage groups will be merged from the set of pseudo-testcross 
#' linkage groups created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG-comb"}, then linkage groups from the set of combined linkage groups 
#' created from the \code{\link{$addBIsnps}} function will be merged.
#' When \code{where = NULL}, linkage groups will be merged in the set of combined linkage groups if available,
#' otherwise they will be merged in the set of pseudo-testcross linkage groups.
#' 
#' In the case when at least one MI linkage group is to be merged with at least one PI linkage group, 
#' the second argument (\code{mergeTo}) specifies whether the merged linkage groups are to be considered
#' as MI linkage groups (\code{maternal}) or PI linkage groups (\code{paternal}). Otherwise,
#' the user will be prompted to input which parental line to merge the linkage groups to.
#' 
#' @usage
#' FSobj$mergeLG(LG, where = NULL, mergeTo = NULL)
#' 
#' @param LG An integer vector specifying the indices of the linkage groups to be merged. 
#' @param where Character vector specifying which list of linkage groups to merge linkage groups from. \code{"LG-pts"} is for 
#' the maternal and paternal linkage groups and \code{"LG-comb"} is for the combined linkage group list.
#' @param mergeTo Character value specifying whether the merged linkage groups are to be considered
#' as MI linkage groups (\code{maternal}) or PI linkage groups (\code{paternal}). Only applies in the case of 
#' merging pseudo-testcross linkage groups.
#' @name $mergeLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## simulate sequencing data
#' set.seed(8932)
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=TRUE), simplify = FALSE))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' 
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' ## merge linkage groups 1 and 2
#' F1data$mergeLG(LG = c(1,2))
#' 
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' ## merge linkage groups 1 and 3
#' F1data$mergeLG(LG = c(1,3), mergeTo="paternal")
NULL