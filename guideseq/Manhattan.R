
args <- commandArgs(trailingOnly=TRUE)

on_target=args[1] # target seq , not used any more
# input_file = "CRL3128_MKSR_Digenomeseq_buffer_rep1_identified_matched.txt"
# output_file = "output.pdf"
input_file = args[2]
output_file = args[3]
# specifically for guide-seq output
library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggforce)
library(ggrepel)
library(viridis)
message(input_file)
CHANGEseq_matched = read_tsv(input_file,col_names=T)

names(CHANGEseq_matched)[names(CHANGEseq_matched) == "BED_Max.Position"] <- "end"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "BED_Min.Position"] <- "start"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "#BED_Chromosome"] <- "chr"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "bi.sum.mi"] <- "reads"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Cell"] <- "sample"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Site_SubstitutionsOnly.NumSubstitutions"] <- "distance"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Site_SubstitutionsOnly.Sequence"] <- "name"

# CHANGEseq_matched <- CHANGEseq_matched %>%
#     mutate(distance = if_else(is.na(distance),6, distance))
# Load hg38 genome object
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# Get chromosome lengths into hg38_tbl
# Note as.numeric() required to store numbers greater than 2E9
# Filter out chromosome strings that contain '_'
# Make a cumulative sum to get y-coordinate of Manhattan plot
# Select all row except for chr_len (temporary variable)
hg38_tbl = tibble( chr = names(seqlengths(hg38)), chr_len = as.numeric(seqlengths(hg38)) ) %>%
  filter(!str_detect(chr,'_')) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len)


# Make an annotated dataset called circleseq_matched_manhattan_annotated
# Join circleseq_matched table with hg38_tbl made above; Left join means that all rows in left table are joined with something in right
# Arrange by sample, chromosome, start
# Add a variable called BPcum that adds the position to the total calculated above
# Select only sample, chr, start, tot, BPcum, and reads for final table

  circleseq_matched_manhattan_annotated = CHANGEseq_matched %>% 
    filter(chr!="chrM") %>%
    left_join(hg38_tbl, by="chr") %>% 
    mutate(sample = str_replace(sample, "CRL[0-9]{3}_", "")) %>%
    arrange(sample, chr, start) %>%
    mutate(BPcum=start+tot) %>%
    select(sample, name,chr, start, tot, BPcum, reads, distance) 
tryCatch({
  circleseq_matched_manhattan_annotated = circleseq_matched_manhattan_annotated %>% 
    group_by(sample) %>%
    mutate(label_x=ifelse(distance==0, BPcum, NA), label_y=ifelse(distance==0, reads, NA), 
           label_distance=ifelse(distance==0, 0.1*max(reads), NA))

},error=function(cond){
  return (NA)
},warning=function(cond){
  return (NA)
},
,finally={message("Couldn't find on-target sequence")}
)
# when match on-target sequence
    # mutate(label_x=ifelse(name==on_target, BPcum, NA), label_y=ifelse(name==on_target, reads, NA), 
           # label_distance=ifelse(name==on_target, 0.1*max(reads), NA))


# Calculate positions to put labels for X-axis  
#axisdf = circleseq_matched_manhattan_annotated %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) %>% arrange(center)
axisdf = hg38_tbl %>% 
  mutate(center=(tot + lead(tot))/2) %>%
  arrange(center)


# Calculate number of pages
npages = ceiling(length(unique(circleseq_matched_manhattan_annotated$sample))/10)

circleseq_matched_manhattan_annotated %>%
  ggplot(aes(x=BPcum, y=reads, color=chr)) +
  # Show all points
  geom_segment( aes(x=BPcum, y=0, xend=BPcum, yend=reads, color=as.factor(chr)), alpha=0.9, size=0.5) +
  
  #geom_text_repel( aes(x=BPcum, y=reads, label=label_text)) + 
  geom_segment(aes(x=label_x, y=reads+2*label_distance, xend=label_x, yend=reads+0.5*label_distance, color="darkred"), 
               linejoin="mitre", arrow=arrow(length=unit(0.2, "cm")), size=0.5) +
  
    # custom X axis:
  scale_x_continuous( name="Chromosome", label = axisdf$chr, breaks= axisdf$center, limits=c(0, 3188269832) ) +
  scale_y_continuous( name="GUIDE-seq V2 read count", expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Customize the theme:
  theme_classic(base_size=15) +
  theme( legend.position="none",
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         strip.background = element_blank(),
         axis.text.x  = element_text(angle=60, hjust=1.1, vjust=1.1, size=8))
  
ggsave(output_file, width=11, height=4)

