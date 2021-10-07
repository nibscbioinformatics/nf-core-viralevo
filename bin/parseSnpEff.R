require(VariantAnnotation)

parseSnpEff <- function(vcf, ...){

  snpeff_description <- info(header(vcf))["ANN","Description"]
  snpeff_string <- sub("^.+: \\'(.+)\\'$", "\\1", snpeff_description)
  snpeff_fields <- gsub(" ", "",
                        unlist(strsplit(snpeff_string, "\\|")))

  .parseAnno <- function(anno){
    fields <- unlist(strsplit(anno, "\\|"))
    if(length(fields)!=length(snpeff_fields)){
      if(length(fields) == length(snpeff_fields)-1){
        ## when string split if last fied is empty the vector will contain -1
        fields=c(fields,"")
      }
      else {
        stop("the number of fields in the annotation string does not correspond to those declared in the header of the VCF file!")
      }
    }
    return(fields)
  }

  .extractAnno <- function(annolist){
    # take the first one which is according to SNPEff the most severe
    annotation <- annolist[1]
    parsed <- .parseAnno(annotation)
    names(parsed) <- snpeff_fields
    return(parsed)
  }

  newMetaData <- do.call(rbind,
                         lapply(info(vcf)$ANN, .extractAnno))

  elementMetadata(vcf) <- cbind(elementMetadata(vcf),
                                newMetaData)
  return(vcf)
}
