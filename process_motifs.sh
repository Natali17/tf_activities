set -euo pipefail

source ./pipeline_core.sh

# Check if the hg38 genome file is already present.
# Returns 0 (success/True) if the file exists, 1 (failure/False) otherwise.
check_genome() {
    if [[ -f source_data/genome/hg38.fa ]]; then
        echo "Genome (hg38.fa) already exists. Skipping download."
        return 0
    fi
    return 1
}

# Check if the HOCOMOCO motif data has been downloaded and extracted.
check_motifs() {
    if [[ -d source_data/motifs/pwm/ ]]; then
        echo "HOCOMOCO motif data already exists. Skipping download."
        return 0
    fi
    return 1
}

# Check if the SARUS JAR file is present and linked.
check_sarus() {
    if [[ -f app/sarus.jar ]]; then
        echo "SARUS jar already exists. Skipping download."
        return 0
    fi
    return 1
}

stage_00_1() {
  echo "Stage 00.1: Downloading hg38 genome"
  rm -rf source_data/genome
  mkdir -p source_data/genome
  pushd source_data/genome
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      gzip -d hg38.fa.gz
  popd
  echo "Done: hg38 genome downloaded and decompressed"
}

stage_00_2() {
  echo "Stage 00.2: Downloading HOCOMOCO motif data"
  rm -rf source_data/motifs/
  mkdir -p source_data/motifs/
  pushd source_data/motifs
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_annotation.jsonl
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pcm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pfm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_thresholds.tar.gz
      for FN in $(ls *.tar.gz); do tar -zxf "${FN}"; done
  popd
  echo "Done: HOCOMOCO motif data downloaded and extracted"
}

stage_00_3() {
  echo "Stage 00.3: Downloading SARUS jar"
  mkdir -p app
  pushd app
      wget https://raw.githubusercontent.com/autosome-ru/sarus/master/releases/sarus-2.1.0.jar
      ln -s sarus-2.1.0.jar sarus.jar
  popd
  echo "Done: SARUS jar ready"
}

stage_00() {
  echo "Stage 00: Downloading all required data"
  # Only execute stages_00_{1:4} if {check_genome, check_motifs, check_sarus} accordingly fails (i.e., data is missing)
  if ! check_genome; then
      stage_00_1
  fi
  
  if ! check_motifs; then
      stage_00_2
  fi
  
  if ! check_sarus; then
      stage_00_3
  fi
  echo "Stage 00 completed"
}

stage_01() {
  echo "Stage 01: Setting up promoters_tss.bed"
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
  TSS_BED="${1:-${SCRIPT_DIR}/new_tss.bed}"

  mkdir -p ./stages/stage_01.preprocess/
  mkdir -p ./stages/stage_01/

  if [[ ! -f "${TSS_BED}" ]]; then
        echo "Error: ${TSS_BED} not found. Please run annotation.py first."
        exit 1
  fi

  rm -f ./stages/stage_01/promoters_tss.bed
  ln -sf "${TSS_BED}" ./stages/stage_01/promoters_tss.bed

  echo "Done: promoters_tss.bed linked to stage_01"
}

stage_02() {
  echo "Stage 02: Generating flanking regions"
  mkdir -p ./stages/stage_02/
  make_flanks 250u 10d  ./stages/stage_01/promoters_tss.bed  ./stages/stage_02/promoters
  #make_flanks 0    50d  ./stages/stage_01/promoters_tss.bed  ./stages/stage_02/promoters
  echo "Done: Flanking regions created"
}

stage_03() {
  echo "Stage 03: Generating FASTA sequences for flanks"
  mkdir -p ./stages/stage_03/
  flanks_fasta 250u 10d  ./stages/stage_02/promoters  ./stages/stage_03/promoters
  #flanks_fasta 0    50d  ./stages/stage_02/promoters  ./stages/stage_03/promoters
  echo "Done: FASTA sequences ready"
}

stage_04() {
  NUM_THREADS=${1:-10} # 10 threads by default
  echo "Stage 04: Computing motif occupancies using ${NUM_THREADS} threads"

  mkdir -p ./stages/stage_04/
  (
    #motif_occupancies_flanks_cmd 250u 10d  ./stages/stage_03/promoters  promoters
    #motif_occupancies_flanks_cmd 0    50d  ./stages/stage_03/promoters  promoters

    motif_besthits_flanks_cmd 250u 10d  ./stages/stage_03/promoters  promoters
    #motif_besthits_flanks_cmd 0    50d  ./stages/stage_03/promoters  promoters
  ) | parallel -j ${NUM_THREADS}
  echo "Done: Motif occupancies computed"
}

stage_00 # downloading data
stage_01
stage_02
stage_03
stage_04 4 # adjust for number of threads available
