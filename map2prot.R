# +-----------------------------------------------------------------------------------+
#  map2prot v1.0.0
# +-----------------------------------------------------------------------------------+

# Notes:
# 1.  This script maps a set of peptides onto a protein sequence. It generates a chart
#      in which the peptides appear as coloured blocks underneath the protein 
#      sequence. Sequence coverage (%) and various metadata are printed on the chart.
# 2.  The script requires:
#      a. R for Windows;
#      b. package 'seqinR';
#      c. a text file containing the peptide sequences;
#      d. a *.fasta file containing one or more protein sequences, plus the index 
#          number of the sequence you want to use (1 for the first or only sequence, 
#          2 for the second sequence, and so on).
# 3.  If the file containing the peptide sequences is called 'peptides.txt', it will 
#      be interpreted as the MaxQuant output file of the same name. If the file is 
#      called something else, it will be interpreted as a plain list (as would be 
#      uploaded for a UniProt peptide search, for example).
# 4.  Please follow these instructions to run the script:
#      a. Open R.
#      b. If necessary, install package 'seqinR'. To do this, type 
#          'install.packages("seqinr")' into the R console (without the single 
#          quotation marks) and press enter.
#      c. Open the script (File > Open script).
#      d. Scroll down to the section headed '>>> Enter parameters here >>>'.
#      e. Complete the path to the file containing your peptide sequences (the 
#          'peptides.file.path').
#      f. Complete the path to your *.fasta file (the 'fasta.file.path').
#      g. Enter the index number of the protein sequence you want to use (the 
#          'fasta.element').
#      h. If using a MaxQuant 'peptides.txt' file, enter two additional pieces of 
#          information:
#           - The name of an 'experiment' in which the peptides were detected 
#              (intensity > 0). This will have originated from the 'Experiment' 
#              column in the raw data table in MaxQuant. If you used the 'No 
#              fractions' button to fill in the table, then your experiments will have
#              the same names as the respective raw data files.
#           - The identifier of a protein of interest, exactly as it appears in 
#              'peptides.txt' (the 'protein.from.which.peptides.derive'). This could 
#              be a UniProt ID, for example.
#      i. Save a copy of the script if required (File > Save as).
#      j. Run the script (Edit > Run all).
#      k. Review the text output in the R console, checking for any warnings.
#      l. Review the graphical output.
# 5.  Blocks (peptides) are placed on the chart in order of decreasing sequence 
#      length. The rule is that blocks cannot overlap or be contiguous on the same 
#      line. The script will keep trying the next line down until it finds an 
#      available space (maximum of eight lines). If it runs out of lines, it will 
#      return an error message stating that a particular peptide will not appear on 
#      the chart.
# 6.  By default, the script is set up to map a set of bovine serum albumin (BSA) 
#      peptides onto the the sequence of mature BSA (A190T variant). The default input
#      files are from PRIDE project PXD013040, and copies are available via 
#      https://github.com/preston-gw/map2prot/. The BSA A190T sequence originates
#      from Protein Data Bank accession 4F5S.
# 7.  map2prot version 1.0.0 was developed in R version 4.2.1, using seqinR version 
#      4.2-16.
# 8.  To cite this script, please consult 
#      https://github.com/preston-gw/map2prot/blob/main/CITATION.cff.

# Start time
start.time <- as.POSIXlt(Sys.time())
print(start.time) # print

# Options
options(warnPartialMatchAttr = TRUE) # set
options(warnPartialMatchDollar = TRUE) # set
getOption("stringsAsFactors") # confirm

# +--------------------------------+
# | >>> Enter parameters below >>> |
# +--------------------------------+

script <- "map2prot_v1.0.0" # will be printed on chart
peptides.file.path <- 
	"C:\\map2prot\\control_replicate1\\txt\\peptides.txt"
fasta.file.path <- "C:\\map2prot\\4f5s_A.fasta"
fasta.element <- 1 # an integer
experiment <- "072617_GPOSU_03E_DDA1.raw" # ignored for plain lists 
protein.from.which.peptides.derive <- "4F5S" # ignored for plain lists

# +--------------------------------+
# | >>>| End of parameters section |
# +--------------------------------+

# Extract file name from 'peptides.file.path'
bits1 <- unlist(strsplit(x = peptides.file.path, split = "\\", fixed = TRUE))
bits2 <- unlist(strsplit(x = bits1[length(bits1)], split = "/", fixed = TRUE))
peptides.file.name <- bits2[length(bits2)]
#| Both types of file path separator ('\\' or '/') should be tolerated.

# Print a message stating how the file will be interpreted
if(peptides.file.name == "peptides.txt")
print(paste("Your file (", 
		peptides.file.name, 
		") will be interpreted as a MaxQuant output file.", 
		sep = ""), 
	quote = FALSE) else
print(paste("Your file (", 
		peptides.file.name, 
		") will be interpreted as plain list.", 
		sep = ""), 
	quote = FALSE)

# Load packages
library(seqinr)

# Define new rounding function
round2 = function(x, digits) {
	posneg = sign(x)
	z = abs(x) * 10^digits
	z = z + 0.5 + sqrt(.Machine$double.eps)
	z = trunc(z)
	z = z / 10^digits
	z * posneg}
#| Unlike the inbuilt 'round' function, this new function always rounds *up* from a 5.
#|  The code is from https://stackoverflow.com/a/12688836.

# Hour, minutes and seconds for time stamp
if(nchar(as.character(start.time$hour)) == 1)
hour <- paste("0", as.character(start.time$hour), sep = "") else
hour <- as.character(start.time$hour)
if(nchar(as.character(start.time$min)) == 1)
minutes <- paste("0", as.character(start.time$min), sep = "") else
minutes <- as.character(start.time$min)
if(nchar(as.character(floor(start.time$sec))) == 1)
seconds <- paste("0", as.character(floor(start.time$sec)), sep = "") else
seconds <- as.character(floor(start.time$sec))

# Time stamp based on start time
time.stamp <- paste(
	gsub(pattern = "-", 
		replacement = "", 
		x = as.Date(as.POSIXlt(start.time))),
	hour, minutes, seconds, sep = "")
print(time.stamp) # print

# Get peptides via one of two routes
if(peptides.file.name == "peptides.txt") { # 'peptides.txt' route

# Read 'peptides.txt'
peptides <- read.table(file = peptides.file.path, 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Print message
print("--------------------------------------------------------------", quote = FALSE)
print(" Number of records in 'peptides.txt'", quote = FALSE)
print("--------------------------------------------------------------", quote = FALSE)
print(paste("   Before filtering:", nrow(peptides)), quote = FALSE)

# If 'experiment' contains hyphens, convert them to full stops.
experiment <- gsub(x = experiment, pattern = "-", replacement = ".")
#| This is necessary because R does an equivalent conversion when it imports column
#|  headings.

# Filter peptides
peptides <- peptides[grepl(pattern = protein.from.which.peptides.derive, 
		x = peptides$Proteins, 
		fixed = TRUE) & # note to self: is 'fixed = TRUE' really required?
	peptides[, paste("Intensity.", experiment, sep = "")] > 0 &
	(peptides$Reverse != "+" | is.na(peptides$Reverse)), ]

# Print message (continued)
print(paste("   After filtering:", nrow(peptides)), quote = FALSE)
print("     Broken down according to protein group:", quote = FALSE)
print(table(peptides$Proteins))
print("--------------------------------------------------------------", quote = FALSE)

# Prepare character vector
peptides <- as.vector(peptides$Sequence)

} else { # 'plain list' route

# Read plain list and convert to character vector
peptides <- read.table(file = peptides.file.path)
peptides <- as.vector(t(peptides))} # end of 'get peptides' chunk

# Sort peptides into order of decreasing length
peptides <- peptides[order(nchar(peptides), decreasing = TRUE)]

# Load *.fasta file
fasta.object <- read.fasta(fasta.file.path,
	seqtype = "AA",
	as.string = TRUE, # want vector rather than string
	whole.header = TRUE) # default is FALSE

# Extract sequence to which peptides will be mapped
sequence.to.which.peptides.will.be.mapped <- fasta.object[[fasta.element]]

# State how many amino acids
print(paste("Sequence to which peptides will be mapped comprises", 
		nchar(sequence.to.which.peptides.will.be.mapped), 
		"amino acid residues."), 
	quote = FALSE)

# Prepare empty matrices
transparency <- matrix(nrow = nchar(sequence.to.which.peptides.will.be.mapped), 
	ncol = 13)
reservations <- matrix(nrow = nchar(sequence.to.which.peptides.will.be.mapped), 
	ncol = 13)

# Prepare 'warnings' vector
warnings <- numeric(nchar(sequence.to.which.peptides.will.be.mapped))
#| If there isn't enough space in the matrix for a given peptide, this will be 
#|  recording in 'warnings' (0 = no warning, 1 = warning).

# Write data into matrices and 'warnings' vector via sliding window
for(i in 1:length(peptides))
	{for(j in 1:(nchar(sequence.to.which.peptides.will.be.mapped) - 
	(nchar(peptides[i]) - 1)))
	{segment <- substr(sequence.to.which.peptides.will.be.mapped, 
		start = j, 
		stop = j + (nchar(peptides[i]) - 1))

	if(segment == peptides[i])
	
	# Define blocks
	{footprint <- j:(j + (nchar(peptides[i]) - 1))
	oversized.footprint <- ((j - 1):
		(j + length(footprint)))[((j - 1):(j + length(footprint))) %in%
			1:nchar(sequence.to.which.peptides.will.be.mapped)]
	# Allocate an unreserved strip, write data into matrices
	if(length(which(is.na(reservations[oversized.footprint, 8]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 8] <- 1
			transparency[footprint, 8] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 7]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 7] <- 1
			transparency[footprint, 7] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 6]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 6] <- 1
			transparency[footprint, 6] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 5]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 5] <- 1
			transparency[footprint, 5] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 4]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 4] <- 1
			transparency[footprint, 4] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 3]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 3] <- 1
			transparency[footprint, 3] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 2]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 2] <- 1
			transparency[footprint, 2] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 1]))) == 
		length(oversized.footprint))
			{reservations[oversized.footprint, 1] <- 1
			transparency[footprint, 1] <- 1} else
	{print(paste("Warning: not enough space in matrix for peptide number ", 
			i, "!", sep = ""), 
		quote = FALSE)
	print(paste(" [Sequence: ", peptides[i], "]", sep = ""), quote = FALSE)
	print("This peptide will not appear in the chart!", quote = FALSE)
	print("---", quote = FALSE)
	warnings[i] <- 1}}}} # record the fact that a warning was given

# Compute sequence coverage using the 'shadow' cast by the peptides
times.position.mapped <- apply(X = transparency, MARGIN = 1, FUN = sum, na.rm = TRUE)
percent.sequence.coverage <- length(
		times.position.mapped[times.position.mapped > 0]) / 
	nchar(sequence.to.which.peptides.will.be.mapped) * 100
percent.sequence.coverage <- round2(percent.sequence.coverage, digits = 1)
#| In order to compute a reliable sequence coverage, the script requires all relevant
#|  peptides to appear in the 'transparency' matrix (as columns of ones).

# Calculate required number of panels for the chart
panel.count <- ceiling(nchar(sequence.to.which.peptides.will.be.mapped) / 100)

# Choose a graphics device and specify plotting-area dimensions
windows(14, 10)

# Define layout
if(panel.count < 6)
par(mfcol = c(7, 1)) else 
par(mfcol = c(panel.count + 1, 1))

# Specify custom plot margins
par(mar = c(3, 5, 0, 5))

# Prepare the graphic, one panel at a time
for(i in 1:panel.count)
{image(x = 1:nchar(sequence.to.which.peptides.will.be.mapped),
	z = transparency,
	breaks = seq(0, 1, 1), 
	col = rainbow(1, start = 2/3, alpha = 0.5),
	xlim = c(i * 100 - 99.5, 
	i * 100 + 0.5),
	axes = FALSE,
	ann = FALSE,
	add = FALSE) # there is no other image!
axis(side = 1, at = c(
		i * 100, 
		i * 100 - 20, 
		i * 100 - 40, 
		i * 100 - 60, 
		i * 100 - 80, 
		i * 100 - 100),
	cex.axis = 1.2) # extends the axis line
for(j in 1:nchar(sequence.to.which.peptides.will.be.mapped))
	text(x = j, y = 0.8, 
		labels = substr(sequence.to.which.peptides.will.be.mapped, 
		start = j, stop = j),
		cex = 1.2)}

# Add text to chart
mtext(side = 1, text = "Residue number", line = 2.5, cex = 0.8) # axis title
mtext(side = 1, 
	text = paste("Script:", script), # script ID
	line = 5, 
	at = panel.count * 100, # right-hand side
	cex = 0.7,
	adj = 1)
mtext(side = 1, 
	text = time.stamp, # time stamp
	line = 6, 
	at = panel.count * 100, # right-hand side
	cex = 0.7,
	adj = 1)
mtext(side = 1, 
	text = paste("Peptides from:", peptides.file.path), # path to *.txt
	line = 5, 
	at = panel.count * 100 - 100, # left-hand side
	cex = 0.7,
	adj = 0)
mtext(side = 1, 
	text = paste("Map-to sequence:", 
		getName(sequence.to.which.peptides.will.be.mapped)), # sequence name 
	line = 6, 
	at = panel.count * 100 - 100, # left-hand side
	cex = 0.7,
	adj = 0)
mtext(side = 1, 
	text = paste("  Source:", fasta.file.path), # path to *.fasta
	line = 7, 
	at = panel.count * 100 - 100, # left-hand side
	cex = 0.7,
	adj = 0)

# If a 'peptides.txt' file was used, add the experiment name to the chart
if(peptides.file.name == "peptides.txt")
{mtext(side = 1, 
	text = paste("Experiment:", experiment), 
	line = 8, 
	at = panel.count * 100 - 100, # left-hand side
	cex = 0.7,
	adj = 0)}

# If there was enough space in 'transparency' for all relevant peptides, print the 
#  sequence coverage on the chart; otherwise, print a warning.  
if(sum(warnings) == 0)
{mtext(side = 1, 
	text = paste("Sequence coverage: ", percent.sequence.coverage, "%", sep = ""), 
	line = 8, 
	at = panel.count * 100, 
	cex = 0.7,
	adj = 1)} else
{mtext(side = 1, 
	text = "Sequence coverage: cannot compute!", 
	line = 8, 
	at = panel.count * 100, 
	cex = 0.7,
	adj = 1)}

# Print an additional warning in the case that there wasn't enough space for all 
#  relevant peptides
if(sum(warnings) > 0)
{mtext(side = 1, 
	text = "WARNING: one or more peptides were not mapped!", 
	line = 7, 
	at = panel.count * 100 - 50, 
	cex = 0.7,
	adj = 1)}

# Print session info
print(sessionInfo())

# Print system info
print(Sys.info())

# Print run end time
print(Sys.time())

# +-----------------------------------------------------------------------------------+
