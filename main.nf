#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f acer-7a0ffa4.sif ]] ;
          then
            singularity pull acer-7a0ffa4.sif docker://index.docker.io/mpgagebioinformatics/acer:7a0ffa4
        fi

    fi


    if [[ "${params.containers}" == "docker" ]] ; 

      then

        docker pull mpgagebioinformatics/acer:7a0ffa4

    fi

    """

}

process proacer {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    val label
    val paired
    val control
    val treatment
    val paired
    val control_gene
  
  script:
    """
#!/usr/local/bin/Rscript
library(ACER)
setwd("${params.output_acer}")
# first, covert counts.count.txt from mageck to the format that can'be properly analyzed by ACER
# columns should be in the order: sgRNA-gene-initial-depleted-initial-depleted--initial-depleted.. 
# i.e. sgRNA	Gene	Control_Rep1	ToxA_Rep1	Control_Rep2	ToxA_Rep2
count_file = read.table("${params.output_mageck_count}counts.count.txt", header = T, check.names=FALSE)
count_file_c = unlist(strsplit("${control}", ","))
count_file_t = unlist(strsplit("${treatment}", ","))
if (length(count_file_c) != length(count_file_t)) {
    print("contorl and treatment should have the same number of samples")
    } else {
    print("rewrite count table for ACER input: counts.count.acer.${label}.txt")}
    
temp <- cbind(count_file[,count_file_c, drop=FALSE], count_file[,count_file_t, drop=FALSE])          # combine
temp <- temp[, c(matrix(1:ncol(temp), nrow = 2, byrow = T))]                 # then reorder
count_file_acer <- cbind(count_file[,1:2], temp)
write.table(count_file_acer,file = "counts.count.acer.${label}.txt", row.names = F, sep = "\\t", quote = FALSE)
# second, check on if negative contrl gene file will be used or not
if (${params.use_neg_ctrl} == T) {
    negctrlfile = as.data.frame(unlist(strsplit("${control_gene}", ",")))
    colnames(negctrlfile) = "Neg_Control_Genes"
    write.table(negctrlfile, file = "neg.gene.${label}.txt", row.names = F, sep = "\\t", quote = FALSE)
    negctrlfilepath = "${params.output_acer}/neg.gene.${label}.txt"
    } else {
    negctrlfilepath = ""}
#now, as counts.count.acer.${label}.txt is ready, we start ACER
#depending if master library is used or not:
if(${params.using_master_library}==F){
${label}.Data <- DataObj\$new(
                    #masterFiles = , 
                    countFile = "${params.output_acer}/counts.count.acer.${label}.txt",  #make sure from the 3rd column, the columns are altenated as initial-depleted-initial-depleted--initial-depleted.. 
                    negCtrlFile = negctrlfilepath, 
                    #sampleInfoFile = "", #subtype file, we don't have so far
                    hasInitSeq = T)
${label}.Model <- ModelObj\$new(user_DataObj = ${label}.Data,
                    use_neg_ctrl= ${params.use_neg_ctrl},
                    #test_samples="", we don't have subtype
                    use_master_library = F) 
}else{
${label}.Data <- DataObj\$new(
                    masterFiles = ${params.acer_master_library}, 
                    countFile = "${params.output_acer}/counts.count.acer.${label}.txt",  #make sure from the 3rd column, the columns are altenated as initial-depleted-initial-depleted--initial-depleted.. 
                    negCtrlFile = negctrlfilepath, 
                    #sampleInfoFile = "", #subtype file, we don't have so far
                    hasInitSeq = T)
${label}.Model <- ModelObj\$new(user_DataObj = ${label}.Data,
                    use_neg_ctrl= ${params.use_neg_ctrl},
                    #test_samples="", we don't have subtype
                    use_master_library = T) 
}
${label}.Result <- optimizeModelParameters(user_DataObj = ${label}.Data,
                    user_ModelObj = ${label}.Model,
                    ncpus = 8)
writeResObj(${label}.Result)
sessionInfo()
    """
}

workflow images {
  main:
    get_images()
}


workflow {
    if ( ! file("${params.output_acer}").isDirectory() ) {
      file("${params.output_acer}").mkdirs()
    }

    rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
    rows=rows.filter{ ! file( "${params.output_drugz}/${it[0]}.txt" ).exists() }
    label=rows.flatMap { n -> n[0] }
    paired=rows.flatMap { n -> n[1] }
    control=rows.flatMap { n -> n[2] }
    control=control.map{ "$it".replace(".fastq.gz","") }
    treatment=rows.flatMap { n -> n[3] }
    treatment=treatment.map{ "$it".replace(".fastq.gz","") }
    control_gene=rows.flatMap { n -> n[5] }

    proacer( label, paired, control, treatment, paired, control_gene)

}