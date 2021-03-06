% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuse-data.R
\docType{data}
\name{modliststring}
\alias{modliststring}
\title{FUSE model structures (string)}
\format{A data frame with 1248 rows and 9 variables}
\description{
This package contains the list of FUSE model structures implemented by Clark et al. (2011). This is available as numerical codes (\code{modlist}) or strings (\code{modliststring}). The variables are as follows:
\itemize{
  \item mid = model id number
  \item rferr = type of rainfall error
  \item arch1 = architecture of the upper soil layer
  \item arch2 = architecture of the lower soil layer
  \item qsurf = surface runoff
  \item qperc = vertical drainage
  \item esoil = evapotranspiration
  \item qintf = interflow
  \item q_tdh = routing
}
Each row identifies a model structure and each column identifies a different element of the model structure. For each element of the model structure, several model decisions can be made (see Clark et al. 2011 for details):
\itemize{
  \item Rainfall error (rferr) \itemize{
      \item \strong{additive_e} (additive rainfall error) = 11
      \item \strong{multiplc_e} (multiplicative rainfall error) = 12
      }
  \item Architecture of the upper soil layer (arch1) \itemize{
      \item \strong{onestate_1} (Single state: upper layer defined by a single state variable) = 21
      \item \strong{tension1_1} (Separate tension storage: upper layer broken up into tension and free storage) = 22
      \item \strong{tension2_1} (Cascading buckets: tension storage sub-divided into recharge and excess) = 23
      }
  \item Architecture of the lower soil layer (arch2) \itemize{
      \item \strong{fixedsiz_2} (baseflow reservoir of fixed size) = 31
      \item \strong{tens2pll_2} (tension reservoir plus two parallel tanks) = 32
      \item \strong{unlimfrc_2} (baseflow reservoir of unlimited size (0-HUGE), frac rate) = 33
      \item \strong{unlimpow_2} (baseflow reservoir of unlimited size (0-HUGE), power recession) = 34
      %% \item topmdexp_2 (topmodel exponential reservoir (no more in the list) = 35
   }
   \item Runoff (qsurf) \itemize{
       \item \strong{arno_x_vic} (Unsaturated zone Pareto: ARNO/Xzang/VIC parameterization (upper zone control)) = 41
       \item \strong{prms_varnt} (Unsaturated zone linear: PRMS variant (fraction of upper tension storage)) = 42
       \item \strong{tmdl_param} (Saturated zone topographic: TOPMODEL parameterization (only valid for TOPMODEL qb)) = 43
       }
   \item Percolation (qperc) \itemize{
       \item \strong{perc_f2sat} (Drainage above field capacity: water from (field cap to sat) avail for percolation) = 51
       \item \strong{perc_lower} (Gravity drainage: percolation defined by moisture content in lower layer (SAC)) = 52
       \item \strong{perc_w2sat} (Saturated zone control: water from wilt point to saturation is available for percolation) = 53
       }
  \item Evaporation (esoil) \itemize{
      \item \strong{rootweight} (root weighting) = 61
      \item \strong{sequential} (sequential evaporation model) = 62
   }
   \item Interflows (qintf) \itemize{
       \item \strong{intflwnone} (interflow denied) = 71
       \item \strong{intflwsome} (interflow allowed) = 72
   }
   \item Routing (q_tdh) \itemize{
       \item \strong{no_routing} (routing denied) = 81
       \item \strong{rout_gamma} (routing allowed using Gamma distribution with shape parameter = 2.5) = 82
       }
 }
 The function \code{readmd()} returns the model components corresponding to a given id. A similar function (\code{readmd2var()}) reads the components as either character or numeric vector.
On Clark's model list, the parent models (Topmodel, VIC, PRMS and Sacramento) have id = 60, 230, 342 and 426 respectively.
For example, model 5 is identified by the following strings:
data(modliststring)
modliststring[5,]
}
\details{
Data set: modliststring
}
