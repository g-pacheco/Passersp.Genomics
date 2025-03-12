#!/bin/bash
set -e

### Script to parse PALEOMIX summary files into a TSV file.
# Version: v0.2.2
# Usage: paleomix_summary2tsv.sh [OPTIONS] PALEOMIX_results_folder1 [PALEOMIX_results_folder2] [...] > stats.tsv
#    PALEOMIX_results_folder1 = folder(s) with PALEOMIX results



#################
### Variables ###
#################

# Check if GNU parallel is installed
if [[ ! `type -P parallel` ]]; then
    echo "ERROR: GNU parallel binary not found in \$PATH."
    exit -1
fi

# Default arguments
N_THREADS=1
MIN_COV=3
N_START=10
K=100
SAMPLES_TO_INCLUDE=''
SUBSETS_TO_PLOT=''
# Parse command line arguments
SHORT=t:c:n:k:
LONG=n_threads:,min_coverage:,n_start:,n_clusters:,samples:,heatmap:
TMP=`mktemp /tmp/paleomix_summary2tsv.$USER.XXXXX`

# -temporarily store output to be able to check for errors
# -activate advanced mode getopt quoting e.g. via “--options”
# -pass arguments only via   -- "$@"   to separate them correctly
PARSED=$(getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    # e.g. $? == 1
    #  then getopt has complained about wrong arguments to stdout
    exit -1
fi
# use eval with "$PARSED" to properly handle the quoting
eval set -- "$PARSED"

# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        -t|--n_threads)
            N_THREADS="$2"
            shift 2
            ;;
        -c|--min_coverage)
            MIN_COV="$2"
            shift 2
            ;;
#        -n|--n_start)
#            N_START="$2"
#            shift 2
#            ;;
#	-k|--n_clusters)
#            K="$2"
#            shift 2
#            ;;
        --samples)
            SAMPLES_TO_INCLUDE="$2"
            shift 2
            ;;
	--heatmap)
	    SUBSETS_TO_PLOT="$2"
	    shift 2
	    ;;
        --)
            shift
            break
	    ;;
        *)
            echo "ERROR: wrong option provided ($1)."
            exit -1
            ;;
    esac
done

case $N_THREADS in
    ''|*[!0-9]*) echo -e "ERROR: invalid number of threads. Check your command:\n\tpaleomix_summary2tsv.sh N_THREADS PALEOMIX_results_folder1 [PALEOMIX_results_folder2] [...] > stats.tsv" >&2; exit -1 ;;
    0*) echo -e "ERROR: number of threads must be > 0." >&2; exit -1 ;;
    *) ;;
esac

# Parse input files to include
find $@ -maxdepth 1 -name '*.*' | 
if [[ -s $SAMPLES_TO_INCLUDE ]]; then
    fgrep -wf $SAMPLES_TO_INCLUDE
else
    cat
fi > $TMP.files

# If defined, exclude all files not in SUBSETS_TO_PLOT
if [[ $SUBSETS_TO_PLOT ]]; then
    # Get SUBSETS present
    grep -E '*.coverage' $TMP.files | xargs -I {} basename {} '.coverage' | sort | awk 'BEGIN{id="#"} { if(gsub(id"[.]","",$1)==0) id=$1; else print $1}' | sort -u > $TMP.subsets
    # Remove SUBSETS
    echo $SUBSETS_TO_PLOT | tr "," "\n" > $TMP.tmp
    fgrep -vwf $TMP.tmp $TMP.subsets > $TMP.subsets_rem
    fgrep -vwf $TMP.subsets_rem $TMP.files > $TMP.tmp
    mv $TMP.tmp $TMP.files
fi



#################
### Functions ###
#################
# print the header (the first line of input)
# and then run the specified command on the body (the rest of the input)
# use it in a pipeline, e.g. ps | body grep somepattern
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
export -f body


########################################
### Convert *.summary files into TSV ###
########################################
awk '
    # Lines to skip
    /^#/ || $1=="Target" || $2!="*" || $4~/_(trash|collapsed)/ {
        next;
    }

    # If first file in list, store header
    NR==FNR {
        header = header"\t"$4;
    }

    # If first record, print Sample ID
    $4=="lib_type" {
        printf "\n"$1;
    }

    # Print value
    {
        printf "\t"$5;
    }

    # Print header
    END {
        print "\nSample_ID"header;
    }
' `grep -E '*.summary' $TMP.files` | tac | awk 'NF' | body sort -k 1,1V



# If no N_START, do not run clustering/heatmap
if [[ $N_START -eq 0 ]]; then
    exit 0
fi



########################
### Find BED regions ###
########################
grep -E '*.coverage' $TMP.files | xargs -I {} basename {} '.coverage' | sort | awk 'BEGIN{id="#"} { if(gsub(id"[.]","",$1)==0) id=$1; else print $1}' | sort -u > $TMP.subsets
N_JOBS=$((`awk 'END{print NR}' $TMP.subsets`))

# Check if number of threads is appropriate
if [[ $N_JOBS -eq 0 ]]; then
    echo "INFO: no SUBSETS found; plotting only overall coverage." >&2
    N_JOBS=1
    MAX_THREADS=1
else
    echo -e "INFO: $N_JOBS SUBSET(s) found: "`perl -n -e 'chomp; push(@s, $_); END{print join(", ",@s)}' $TMP.subsets` >&2
    MAX_THREADS=$((N_JOBS*N_START))
fi

if [[ $N_THREADS -gt $MAX_THREADS ]]; then
    echo "ERROR: too many threads specified! Only $MAX_THREADS can actually be used at the most..." >&2
    exit -2
fi



######################
### Make TSV files ###
######################
# Function to parse PALEOMIX coverage files to TSV
parse_coverage() {
    awk '/#/ || /Coverage/ || $2!="*" || $4=="*" {cnt=1; next} NR==FNR{if(cnt==1)h=$4; else h=h"\t"$4} cnt++==1{printf "\n"$1} {printf "\t"$14} END{print "\n"h}' `grep -E "*.$1.coverage" $2` | awk 'NF' | tac | LC_COLLATE=C body sort -k 1,1V
}
# Function to, given a threshold, convert values to presence/absence
# Usage: qual2quant COV FILE
qual2quant(){
    cat $2 | awk -v min_cov=$1 'NR==1{print; next} { printf $1; for(i=2;i<=NF;i++) printf "\t"($i >= min_cov ? 1 : 0); print "" }'
}
export -f parse_coverage qual2quant

# Convert depth from PALEOMIX format to TSV
# Change values to percentage [0,100] (in PALEOMIX is [0,1])
# ALL
awk '/#/{next} $2!="*" || $4!="*"{cnt=1; id=FILENAME; gsub(".*/","",id); gsub(/\.depths/,"",id); if(!/MaxDepth/ || NR!=FNR) next; cnt--} cnt++==1{printf "\n"id} {for(i=7;i<=NF;i++) printf (cnt==1 && i==7 ? "" : "\t")(cnt==1 ? $i : $i*100)}' `grep -E '*.depths' $TMP.files | fgrep -vwf $TMP.subsets` | body sort -k 1,1V > ALLc.coverage.tsv
# ALL (non-cumulative)
awk 'BEGIN{OFS="\t"} NR==1{gsub("MD_","D_",$0)} NR>1{ for(i=2;i<NF;i++){$i-=$(i+1)} } {print}' ALLc.coverage.tsv > ALL.coverage.tsv
# SUBSETS (converted to presence/absence; assuming depth >=3)
parallel -j $N_THREADS "parse_coverage {} $TMP.files | qual2quant $MIN_COV - > {}.coverage.tsv" :::: $TMP.subsets



exit 0

####################
### Plot heatmap ###
####################
plot_heatmap() {
cat <<EOF | R --vanilla --slave --args $1 $2 $3 $4 $5
  library(stringi, quiet=TRUE)
  library(bigmemory, quiet=TRUE)
  library(bigpca, quiet=TRUE)
  library(pheatmap, quiet=TRUE)
  args <- commandArgs(trailingOnly = TRUE); cpu=as.numeric(args[1]); id=args[2]; k=as.numeric(args[3]); nstart=as.numeric(args[4]); source_path=args[5]
  cat("\n====> Run:", id, file=stderr(), fill=TRUE)
  # Check number of CPUs
  if(cpu < 1)
    cpu = 1
  # Define unique string
  rnd_id <- stri_rand_strings(n=1, length=20, pattern="[A-Za-z0-9]")
  data <- read.big.matrix(paste(id,"tsv",sep="."), sep="\t", header=TRUE, has.row.names=TRUE, type="short")#, backingfile=paste(rnd_id,"m","bck", sep="."), descriptorfile=paste(rnd_id,"m","dsc", sep="."));
  cat("Number of rows:", nrow(data), file=stderr(), fill=TRUE)
  cat("Number of columns:", ncol(data), file=stderr(), fill=TRUE)
  # Transverse
  data <- big.t(data, name=rnd_id, delete.existing=FALSE);
  cat("Number of transversed rows:", nrow(data), file=stderr(), fill=TRUE)
  cat("Number of transversed columns:", ncol(data), file=stderr(), fill=TRUE)

  # Check if enough data for a heatmap and hierarchical clustering
  if(nrow(data) > 2) {
    # Check if enough data for k-means clustering
    if(nrow(data) > k &&
       id != 'ALL.coverage' &&
       id != 'ALLc.coverage') {
      library(biganalytics)
      source(paste(source_path,"kmeans.R",sep="/"), chdir=TRUE)
      set.seed(12345678);

      uniq <- data[!duplicated(as.matrix(data)), ]
      if(nrow(uniq) < k){
        k = nrow(uniq);
        nstart = 1;
      }

      cpu = min(cpu, nstart)
      if(cpu > 1){
        library(doMC)
        registerDoMC(cpu);
      }

      cat("Number of unique rows:", nrow(uniq), file=stderr(), fill=TRUE)
      cat("Number of Ks:", k, file=stderr(), fill=TRUE)
      cat("Number of independent starts:", nstart, file=stderr(), fill=TRUE)
      cat("Number of threads:", cpu, file=stderr(), fill=TRUE)
      km <- bigkmeans(data, centers=k, iter.max=1000, nstart=nstart)
      # Save clusters info to file
      cat("Saving clustering", file=stderr(), fill=TRUE)
      write.table(matrix(km[['cluster']], ncol=1, dimnames=list(rownames(data),"Cluster")), file=paste(id,"km_cluster","tsv",sep="."), quote=FALSE, sep="\t")

      # Prepare K-means clustering matrix for plotting
      col_names <- colnames(data)
      data <- km[["centers"]];
      clust_sizes <- table(km[["cluster"]]);
      rownames(data) <- sprintf("Cluster: %s Size: %d", names(clust_sizes), clust_sizes);
      colnames(data) <- col_names
      # Treat all negative values as zero
      data[data < 0] <- 0
      # Save raw data to file
      write.table(data, file=paste(id,"km","tsv",sep="."), quote=FALSE, sep="\t")

      # Calculate weighted data (by cluster size)
      data_weighted <- data * sqrt(as.vector(clust_sizes))
      # Save weighted data matrix to file
      write.table(data_weighted, file=paste(id,"km_weighted","tsv",sep="."), quote=FALSE, sep="\t")
      # Get distance
      dist_weighted <- dist(t(data_weighted), method="manhattan")
      # Sort data matrix by size of cluster
      data <- data[rev(order(clust_sizes)),]
    }else{
      cat("WARN: Not enough regions in", id, "for k-means clustering.", file=stderr(), fill=TRUE)
      data <- as.matrix(data)
      dist_weighted <- "manhattan"
    }

    # Plot Heatmap
    nc <- ncol(data)
    nr <- nrow(data)
    pheatmap(data, border_color="black", clustering_distance_cols=dist_weighted, clustering_distance_rows="manhattan", treeheight_row=max(10,floor(nr*0.5)), treeheight_col=max(10,floor(nc*0.5)), cellwidth=10, cellheight=10, cutree_rows=max(5,floor(nr*20/1000)), cutree_cols=max(5,floor(nc*20/1000)), filename=paste(id,"TREE","pdf",sep="."));
    pheatmap(data, border_color="black", cluster_rows=FALSE, cluster_cols=FALSE, cellwidth=10, cellheight=10, filename=paste(id,"pdf",sep="."));
  }else{
    cat("WARN: Not enough regions in", id, "to plot heatmap.", file=stderr(), fill=TRUE)
  }

  # Clean-up
  for (file in list.files(pattern=rnd_id)) {
    file.remove(file)
  }
EOF
}
export -f plot_heatmap
\parallel -j $N_JOBS "plot_heatmap $((N_THREADS/N_JOBS)) {.} $K $N_START `dirname $0`" ::: *.coverage.tsv

# Remove files that do not make sense
rm -f ALL.coverage.TREE.pdf ALLc.coverage.TREE.pdf
# Remove temporary files
rm -f $TMP.* Rplots.pdf
