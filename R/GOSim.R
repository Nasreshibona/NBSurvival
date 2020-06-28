# Copyright (C) 2020  Hocine Bendou <hocine@sanbi.ac.za>
#                     Abdulazeez Giwa <3901476@myuwc.ac.za>
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>

library(AnnotationHub)
library(GOSemSim)

# Build annotation data. We using Gene Symbol as keytype (Default is Gene entriz)
hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)

# Upregulated genes. Some of the genes don't have GO information (They are removed from the list)
genes <- c("EVX2","NHLH2","POU6F2","HOXD10","RTL1","MAPK15","LGR5","PRSS12","STRA6")
mgeneSim(genes, semData=hsGO,measure="Wang",combine="BMA", verbose=FALSE)

# OUTPUT
#         EVX2 NHLH2 POU6F2 HOXD10 MAPK15  LGR5 PRSS12 STRA6
# EVX2   1.000 1.000  1.000  1.000  0.214 0.214  0.119 0.139
# NHLH2  1.000 1.000  0.643  1.000  0.484 0.581  0.145 0.169
# POU6F2 1.000 0.643  1.000  1.000  0.197 0.224  0.145 0.169
# HOXD10 1.000 1.000  1.000  1.000  0.214 0.214  0.119 0.139
# MAPK15 0.214 0.484  0.197  0.214  1.000 0.483  0.226 0.169
# LGR5   0.214 0.581  0.224  0.214  0.483 1.000  0.145 0.169
# PRSS12 0.119 0.145  0.145  0.119  0.226 0.145  1.000 0.100
# STRA6  0.139 0.169  0.169  0.139  0.169 0.169  0.100 1.000

# Downregulated genes
genes <- c("CYP17A1","MYH14","LRRTM3","GRIN3A","HS3ST5","NBAS","FNDC9","HIST1H1E","CRYAB","NXPH3","MYL3","CMYA5","AMIGO2","SIK1B","EDIL3","UBC")
mgeneSim(genes, semData=hsGO,measure="Wang",combine="BMA", verbose=FALSE)

# OUTPUT
#          CYP17A1 MYH14 GRIN3A HS3ST5  NBAS HIST1H1E CRYAB NXPH3  MYL3 CMYA5 AMIGO2 EDIL3   UBC
# CYP17A1    1.000 0.199  0.333  0.332 0.481    0.425 0.417 0.477 0.253 0.590  0.590 0.262 0.427
# MYH14      0.199 1.000  0.269  0.251 0.419    0.287 0.349 0.234 0.331 0.508  0.508 0.133 0.301
# GRIN3A     0.333 0.269  1.000  0.520 0.815    0.543 0.677 0.477 0.370 1.000  1.000 0.262 0.595
# HS3ST5     0.332 0.251  0.520  1.000 0.815    0.506 0.616 0.477 0.367 1.000  1.000 0.262 0.575
# NBAS       0.481 0.419  0.815  0.815 1.000    0.815 0.815 0.383 0.482 0.815  0.815 0.214 0.815
# HIST1H1E   0.425 0.287  0.543  0.506 0.815    1.000 0.674 0.477 0.388 1.000  1.000 0.262 0.868
# CRYAB      0.417 0.349  0.677  0.616 0.815    0.674 1.000 0.477 0.428 1.000  1.000 0.262 0.699
# NXPH3      0.477 0.234  0.477  0.477 0.383    0.477 0.477 1.000 0.477 0.477  0.477 0.477 0.477
# MYL3       0.253 0.331  0.370  0.367 0.482    0.388 0.428 0.477 1.000 0.602  0.602 0.590 0.427
# CMYA5      0.590 0.508  1.000  1.000 0.815    1.000 1.000 0.477 0.602 1.000  1.000 0.262 1.000
# AMIGO2     0.590 0.508  1.000  1.000 0.815    1.000 1.000 0.477 0.602 1.000  1.000 0.262 1.000
# EDIL3      0.262 0.133  0.262  0.262 0.214    0.262 0.262 0.477 0.590 0.262  0.262 1.000 0.262
# UBC        0.427 0.301  0.595  0.575 0.815    0.868 0.699 0.477 0.427 1.000  1.000 0.262 1.000

# NHLH2 with AMIGO2, CMYA5, CRYAB and NBAS. NHLH2 is up and the other genes are down.
genes <- c("NHLH2","AMIGO2","CMYA5","CRYAB","NBAS")
mgeneSim(genes, semData=hsGO,measure="Wang",combine="BMA", verbose=FALSE)

# OUTPUT
#        NHLH2 AMIGO2 CMYA5 CRYAB  NBAS
# NHLH2  1.000  1.000 1.000 0.664 0.815
# AMIGO2 1.000  1.000 1.000 1.000 0.815
# CMYA5  1.000  1.000 1.000 1.000 0.815
# CRYAB  0.664  1.000 1.000 1.000 0.815
# NBAS   0.815  0.815 0.815 0.815 1.000

genes <- c("CYP17A1","MYH14","LRRTM3","GRIN3A","HS3ST5","NBAS","FNDC9","HIST1H1E","CRYAB",
           "NXPH3","MYL3","CMYA5","AMIGO2","SIK1B","EDIL3","UBC",
           "EVX2","NHLH2","POU6F2","HOXD10","RTL1","MAPK15","LGR5","PRSS12","STRA6")
mgeneSim(genes, semData=hsGO,measure="Wang",combine="BMA", verbose=FALSE)
