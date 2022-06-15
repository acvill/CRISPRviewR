read_minced <- function(dir) {

  seqs <- read_tsv(file = paste0(dir,"/minced_out.txt"),
                   col_names = "mess", comment = "--") %>%
    filter(!grepl('POSITION', mess)) %>%
    filter(!grepl('Repeats', mess)) %>%
    mutate(mess = str_replace(mess, "CRISPR ", "CRISPR")) %>%
    mutate(mess = str_replace(mess, "Range.*", "")) %>%
    group_by(grp = cumsum(str_detect(mess, "CRISPR"))) %>%
    mutate(array = dplyr::first(mess)) %>%
    filter(!grepl('CRISPR', mess)) %>%
    ungroup %>%
    separate(col = mess, remove = T, sep = "\t",
             into = c("start","drop1","rep","spacer","drop2")) %>%
    mutate(array = str_trim(array)) %>%
    mutate(sample = as.character(dir)) %>%
    dplyr::select(sample, array, start, rep, spacer) %>%
    mutate(start = as.character(start))

  coords <- read_tsv(file = paste0(dir,"/minced_out.gff"),
                     comment = "##",
                     col_names = c("contig","drop1",
                                   "drop2","start","end","drop3",
                                   "drop4","drop5","parse")) %>%
    dplyr::select(-contains("drop")) %>%
    filter(!grepl("rpt_type", parse)) %>%
    separate(col = parse, sep = ";", into = c("array","id"), remove = TRUE) %>%
    mutate(array = str_remove(pattern = "Parent=", string = array)) %>%
    mutate(sample = as.character(dir)) %>%
    dplyr::select(sample, array, contig, start, end) %>%
    mutate(start = as.character(start))

  left_join(seqs, coords, by = c("sample", "array", "start")) %>%
    mutate(start = as.numeric(start)) %>%
    dplyr::select(sample, array, rep, start, end, spacer, contig)

}
