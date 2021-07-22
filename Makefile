# The complete makefile will run to the compile step
all: compile

# Remove previous processing output
clean:
	rm -r output

# Identify script directory
source = code

# Demultiplex reads
demux.out = output/demux
demux: $(demux.out)/
$(demux.out)/: $(source)/canopy.pheniqs.config.json
	mkdir -p $(demux.out)
	pheniqs mux -R $(demux.out)/demux.report.txt -c $< --base-input data/raw --base-output $(demux.out)

# Trim primers / adapters with cutadapt and SeqPurge
trim.out = output/trim
trim: demux $(trim.out)/
$(trim.out)/: $(source)/trim.sh
	$<

# Denoise with dada2
denoise.out = output/denoise/seq.tab.rds
denoise: trim $(denoise.out)
$(denoise.out): $(source)/denoise.R
	Rscript $<

# Trims with ITSx, removes chimeras, postcluster with LULU, filter with PERFect,
# adds taxonomy, and saves output as a list of phyloseq objects
compile.out = output/compile/phy.list.rds
compile: denoise $(compile.out)
$(compile.out): $(source)/compile.R
	Rscript $<
