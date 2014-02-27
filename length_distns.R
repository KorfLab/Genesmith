#============================#
# <length_distns.R>          #
# -  Plots gene length distn #
# -  Plots CDS/Intron distn  #
#============================#

# Import Concatenated Output files from <kogs_distns.pl>
KOG_DATA     = read.table("~/Work/KorfProgramProjects/Genesmith/bin/kog_gene_feat.txt")
CDS_DATA     = read.table("~/Work/KorfProgramProjects/Genesmith/bin/kog_cds_len.txt")
IN_DATA      = read.table("~/Work/KorfProgramProjects/Genesmith/bin/kog_in_len.txt")

# Name Columns
names(KOG_DATA) = c("taxa", "dna_length", "aa_length", "ex_quant", "in_quant")
names(CDS_DATA) = c("taxa", "rank", "length")
names(IN_DATA)  = c("taxa", "rank", "length")

# Split up Datasets Based on Species
kog_genomic  = split(KOG_DATA$dna_length, KOG_DATA$taxa)
kog_protein  = split(KOG_DATA$aa_length, KOG_DATA$taxa)
kog_ex_quant = split(KOG_DATA$ex_quant, KOG_DATA$taxa)
kog_in_quant = split(KOG_DATA$in_quant, KOG_DATA$taxa)
kog_cds      = split(CDS_DATA$length, CDS_DATA$taxa)
kog_in       = split(IN_DATA$length, IN_DATA$taxa)

# KOG Taxa/Color Code
kog_colors  = c("green", "red", "blue", "black", "purple", "orange")
kog_species = c("A.thaliana", "C.elegans", "D.melanogaster", "H.sapiens", "S.cerevisiae", "S.pombe") 

#------------------------------#
# KOG Gene Length Distribution #
#------------------------------#
densities  = lapply(kog_genomic, density, na.rm = TRUE)
length(densities) == length(kog_genomic)    # sanity check

table(sapply(densities, class))
xrange    = range(sapply(densities, function(d) range(d$x)))
yrange    = range(sapply(densities, function(d) range(d$y)))
lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Longer Plot
plot(densities[[1]], col = kog_colors[1],
     xlim = lrg_xlim, ylim = yrange,
     main = "KOGs Gene Length Distribution")
# Plot Zoomed in
plot(densities[[1]], col = kog_colors[1],
     xlim = small_xlim, ylim = yrange,
     main = "KOGs Gene Length Distribution")

# Add Other Density Plots
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       densities[-1], kog_colors[-1])

# Add Color Key in the top right
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)


#-----------------------------#
# KOG CDS Length Distribution #
#-----------------------------#
densities  = lapply(kog_cds, density, na.rm = TRUE)
length(densities) == length(kog_cds)    # sanity check

table(sapply(densities, class))
xrange    = range(sapply(densities, function(d) range(d$x)))
yrange    = range(sapply(densities, function(d) range(d$y)))
lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Plot and Label
plot(densities[[1]], col = kog_colors[1],
     xlim = lrg_xlim, ylim = yrange,
     main = "KOGs CDS Length Distribution")
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       densities[-1], kog_colors[-1])
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)

#---------------------------------#
# KOGs CDS Frequency Distribution #
#---------------------------------#
densities  = lapply(kog_ex_quant, density, na.rm = TRUE)
length(densities) == length(kog_ex_quant)    # sanity check

table(sapply(densities, class))
xrange    = range(sapply(densities, function(d) range(d$x)))
yrange    = range(sapply(densities, function(d) range(d$y)))
lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Plot and Label
plot(densities[[1]], col = kog_colors[1],
     xlim = lrg_xlim, ylim = yrange,
     main = "KOGs CDS Frequency Distribution")
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       densities[-1], kog_colors[-1])
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)


#--------------------------------#
# KOG Intron Length Distribution #
#--------------------------------#
densities  = lapply(kog_in, density, na.rm = TRUE)
length(densities) == length(kog_in)    # sanity check

table(sapply(densities, class))
xrange    = range(sapply(densities, function(d) range(d$x)))
yrange    = range(sapply(densities, function(d) range(d$y)))
lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Plot and Label
plot(densities[[1]], col = kog_colors[1],
     xlim = lrg_xlim, ylim = yrange,
     main = "KOGs Intron Length Distribution")
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       densities[-1], kog_colors[-1])
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)


#-----------------------------------#
# KOG Intron Frequency Distribution #
#-----------------------------------#
densities  = lapply(kog_in_quant, density, na.rm = TRUE)
length(densities) == length(kog_in_quant)    # sanity check

table(sapply(densities, class))
xrange    = range(sapply(densities, function(d) range(d$x)))
yrange    = range(sapply(densities, function(d) range(d$y)))
lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Plot and Label
plot(densities[[1]], col = kog_colors[1],
     xlim = lrg_xlim, ylim = yrange,
     main = "KOGs Intron Frequency Distribution")
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       densities[-1], kog_colors[-1])
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)

#---------------------------------#
# KOG Protein Length Distribution #
#---------------------------------#
aa_densities = lapply(kog_protein, density, na.rm = TRUE)
length(aa_densities) == length(kog_protein)    # sanity check

# Set the Range
table(sapply(aa_densities, class))
aa_xrange    = range(sapply(aa_densities, function(d) range(d$x)))
aa_yrange    = range(sapply(aa_densities, function(d) range(d$y)))
aa_lrg_xlim   = c(-12486.27,100000)  # temp xlim with large range
aa_small_xlim = c(-1000,10000)       # temp xlim with smaller range

# Plot and Label
plot(aa_densities[[1]], col = kog_colors[1],
     xlim = aa_lrg_xlim, ylim = aa_yrange,
     main = "KOGs Protein Length Distribution")
mapply(function(d, kog_colors)
  lines(d, col = kog_colors),
       aa_densities[-1], kog_colors[-1])
legend("topright", kog_species, col = kog_colors, lty = 1, cex = .70)

