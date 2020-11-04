#' The function extract_region uses a reference fasta file to set up a custom BLAST database.
#' It then queries all fasta files in the query folder, and returns regions with BLAST matches to the database.
#'
#' Before use, BLAST+ must be installed on the local machine.
#'
#' @param ref Fasta file name: the file containing the fasta-format sequences which will be set up as a database.
#' @param qry_fld Folder name: the query folder, containing all sequences (fasta format) that you wish to extract regions from. Files in the folder should have only the suffixes ".fasta", ".fna", ".fa" or ".fas"
#' @param temp_dir Folder name: the temporary folder in which to store files
#' @param len_thresh Integer: A minimum sequence length to retain for BLAST matches.
#'
#' @return A data frame summarising the BLAST results from the query
extract_region <- function(ref=NULL,qry_fld=NULL,temp_dir=NULL,len_thresh=NULL) {
  if (is.null(ref)) stop("Please provide a fasta file containing reference sequence data")
  if (is.null(temp_dir)) temp_dir=tempdir()

  ref=normalizePath(ref)
  temp_dir=normalizePath(temp_dir)

  print("This function extracts DNA regions that best match the provided reference sequence from all sequences in a specified folder, using BLASTN.")
  print("The folder should contain files with the extensions .fa .fas .fna .fasta")
  print(paste0("Temporary files will be written to this location: ",temp_dir))

  lst=list.files(path=qry_fld,pattern=".fa$|.fas$|.fasta|.fna",full.names = T)
  nm=length(lst)
  print(paste0("The data folder contains ", nm," files to be queried."))

  print("Reading reference fasta file.")
  tryCatch(expr={fas=Biostrings::readDNAStringSet(ref)},error=function(e){stop("ERROR: Please provide a valid fasta-format reference file.")})
  print(paste0("The reference file contains ",length(fas)," sequence(s)."))

  print("Generating BLAST database")
  pt=normalizePath(paste0(temp_dir,"/sequences"),winslash = "\\",mustWork = F)
  print(pt)
  pt2=normalizePath(paste0(temp_dir,"/results"),winslash = "\\",mustWork = F)
  print(pt2)
  system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ",ref," -out ",pt))
  print("Done.")

  print("BLASTing all query sequences, and extracting the best hits")

  res=list()
  for (x in 1:nm){
    print(paste0("BLASTing sequence ",x))
    print(paste0("blastn -num_threads 4 -query ",normalizePath(lst[x])," -db ",pt," -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out ",pt2))

    system(paste0("blastn -num_threads 4 -query ",normalizePath(lst[x])," -db ",pt," -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" -out ",pt2))

    #
    if(file.info(pt2)$size>0){
      tab=read.table(pt2,header=F,stringsAsFactors = F)
      colnames(tab)=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")
      res[[x]]= tab %>% dplyr::group_by(qaccver) %>% dplyr::summarise(best_hit=saccver[order(evalue,-rank(bitscore))[1]],qseq=qseq[order(evalue,-rank(bitscore))[1]],score=evalue[order(evalue,-rank(bitscore))[1]],length=length[order(evalue,-rank(bitscore))[1]],file=basename(lst[x]))
    }else{
      print(paste0("Sequence ",x," (",basename(lst[x]),") has no BLAST hits"))
    }
  }
  res=do.call(rbind,res)
  if(!is.null(len_thresh)){
    res=res[res$length>len_thresh,]
  }
  return(res)
}

#' The function align_df aligns sequences returned by extract_region using MAFFT, and assigns an allele number to each sequence.
#'
#' Before use, MAFFT must be installed on the local machine.
#'
#' @param df dataframe: the output of the function extract_region.
#' @param temp_dir Folder name: the temporary folder in which to store files
#'
#' @return A data frame with aligned sequences and allele numbers
align_df=function(df=NULL,temp_dir=NULL){
  if (is.null(temp_dir)) temp_dir=tempdir()

  print("This function aligns the sequences in the output from function extract_region using system calls to mafft.")
  print(paste0("Temporary files will be written to this location: ",temp_dir))
  print("Writing sequence data to disk")
  seqinr::write.fasta(sequences = as.list(gsub("-| ","",df$qseq)),names = as.list(df$file),file.out = paste0(temp_dir,"/aln.fasta"),open="w")
  print("Calling mafft")
  setwd(temp_dir)
  print(paste0("mafft --adjustdirection --thread 4 aln.fasta > aln2.fasta"))
  system(paste0("mafft --adjustdirection --thread 4 aln.fasta > aln2.fasta"),intern=T)
  dna=Biostrings::readDNAStringSet("aln2.fasta")
  df=data.frame(file=names(dna),id=df$qaccver,aln=paste(dna))
  un=unique(df$aln)
  df$allele=match(df$aln,un)
  return(df)
}


#' The function get_source_data parses eutils to return source metadata for a list of Genbank IDs.
#'
#' @param ids character vector: containing all IDs that you wish to extract metadata for.
#'
#' @return A data frame with all source metadata
get_source_data=function(ids=NULL){

  md=lapply(1:length(ids),function(y){
    print(paste0(y," of ",length(ids)))
    URL=paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",ids[y],"&rettype=gb&retmode=text")
    data <- tryCatch(scan(file = URL, what = "", sep = "\n", quiet = TRUE),error=function(cond){NA},warning=function(cond){NA})
    if(!is.na(data)){
      pos1=which(unlist(sapply(1:length(data),function(x) length(grep("     source",data[x]))==1)))+1
      pos2=which(unlist(sapply(1:length(data),function(x) length(grep("                     ",data[x]))==0)) & sapply(1:length(data),function(x) x>pos1))[1]-1
      tmp=lapply(pos1:pos2,function(x){
        tmp=str_split(str_match(data[x], "^                     /(.*?)\"$")[[2]],"=\"",simplify = T)
        dftmp=data.frame(tmp=tmp[2])
        colnames(dftmp)=tmp[1]
        return(dftmp)

      })
      return(cbind(data.frame(id=ids[y],stringsAsFactors = F),bind_cols(tmp)))
    }else{
      return(data.frame(id=ids[y],stringsAsFactors = F))
    }
  })
  md=bind_rows(md)
  return(md[,!grepl("V[0-9]",colnames(md))])

}
