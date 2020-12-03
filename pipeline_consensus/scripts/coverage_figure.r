

# This script should be made so that it takes a number of arguments
# 1. A glob path to unfiltered+filtered depth .tsv's containing the samtools depth columns preceeded by the sample name
#        Or just the files catted together.
# 2. The path to the bed-file containing the primers wanted for visualization.
# 3. The path to the plot produced

args = commandArgs(trailingOnly=TRUE)
print("These are the args")
print(args)


path_catted_coverage = args[1] #  "~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/output/201123/depth_all.tsv"
path_primer_bed = args[2] # "~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/artic-ncov2019_gh_clone/artic-ncov2019/primer_schemes/SARS-CoV-2/V3/nCoV-2019.primer.bed"
path_insert_bed = args[3] # "~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/artic-ncov2019_gh_clone/artic-ncov2019/primer_schemes/SARS-CoV-2/V3/nCoV-2019.insert.bed"
path_out = args[4] # plots/cov.png
cowplot_source_path = args[5]  

if (length(args) > 4) {
    print("assuming cowplot source is given as 5th argument")
    if (!require("cowplot")) {
        install.packages(cowplot_source_path, lib = tail(.libPaths(), 1), repos = NULL, type="source")
    }  
    

}






library(tidyverse)
library(cowplot)



# concatenate and read all depth-files

#cmd = "cat ~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/output/201123/*/*depth.tsv > ~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/output/201123/depth_all.tsv"
#system(cmd)


#coverage = read_tsv("~/GenomeDK/ClinicalMicrobio/faststorage/artic/ivar/output/201123/depth_all.tsv", col_names = c("sample", "filter", "target","position", "coverage"))
coverage = read_tsv(path_catted_coverage, col_names = c("sample", "filter", "target","position", "coverage"))

str(coverage)

coverage$sample %>% table


# read bed-file
bed_primer = read_tsv(path_primer_bed,
                      col_names = c("chrom", "chrom_start", "chrom_end", "name", "score", "strand"))

str(bed_primer)

bed_insert = read_tsv(path_insert_bed,
                      col_names = c("chrom", "chrom_start", "chrom_end", "name", "score", "strand")) %>% 
    rowwise() %>% 
    mutate(y = if_else(score == 1, -1, 1),
           mid = mean(c(chrom_start, chrom_end)))

str(bed_insert)


#bed_list = bed %>% pivot_longer(c(chrom_start, chrom_end), names_repair = "unique")


# plot coverage per sample
a = coverage %>% 
    #filter(sample == "auh-s67_V53490") %>% 
    
    ggplot(aes(position, coverage)) + 
    
   
    #geom_line(data = coverage %>% filter(filter == "before trimming"), aes(position, coverage), linetype = "dashed") + 
    #geom_line(data = coverage %>% filter(filter == "after trimming"), aes(position, coverage), linetype = "solid") + 
    
    # primer horizontal lines
    geom_segment(data = bed_primer, mapping = aes(x = chrom_start, xend = chrom_end, y = 0, yend = 0), color = "gray90", size = 8000, alpha = 0.5) +
    
    
    geom_area(data = coverage %>% filter(filter == "before trimming"), aes(position, coverage), fill = "grey25") + 
    geom_area(data = coverage %>% filter(filter == "after trimming"), aes(position, coverage), fill = "grey") +  
    
    
    
    
    #geom_area(alpha = 0.3) + # time intensive
    #scale_color_manual(values = c("#00bfc4", "#f8766d")) + 
    
    #scale_y_log10() +
    
    facet_grid(sample~.) + 
    theme_minimal() +
    theme(legend.position = "none",
          strip.text.y.right = element_text(angle = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    xlim(0, coverage$position %>% max) +
    labs(); #a


 #ggsave("plots/cov_a.png", height = 10, width = 30, plot = a)    




# blå   00bfc4
# lilla c77cff
# rød   f8766d
# grøn  7cae00
    

# inserts
#b = data = bed_insert %>%
b = bed_insert %>% 
    ggplot() + 
    geom_segment(mapping = aes(x = chrom_start, xend = chrom_end, 
                               y = y, yend = y,
                               color = as.factor(score)), size = 2) + 
    geom_text(aes(x = mid, y = y/4, label = name), size = 3)+
    theme_void() + 
    theme(legend.position = "none") + 
    geom_segment(data = bed_primer, mapping = aes(x = chrom_start, xend = chrom_end, y = 0, yend = 0, color = as.factor(strand)), size = 8000) +
    

    scale_color_manual(values = c("#00bfc4", "#f8766d", "#c77cff", "#7cae00")) +
    
    facet_grid(chrom~.) +
    xlim(0, coverage$position %>% max) + 
    ylim(-1.5,1.5); b



p = plot_grid(b, a, ncol = 1,
              rel_heights = c(1, 10),
              align = "v", axis = "lr")
# p


#ggsave("plots/cov.png", height = 10, width = 30, plot = p)    
ggsave(path_out, height = 10, width = 30, plot = p)    


# deprecated stuff
# b/a +
#     plot_layout(heights = c(1, 20))

#align_plots(b, a, align="h", axis="lr")

#ggdraw(t[[1]]) + draw_plot(t[[2]])

#library(patchwork)
#t = align_patches(a, b)




#
#
#
#
#lst = list()
#for (group in group_split(coverage %>% group_by(sample))) {
#    glimpse(group)
#    
#    sample_name = group$sample %>% unique
#    plot_part = group %>% 
#        #filter(sample == "auh-s67_V53490") %>% 
#        
#        ggplot(aes(position, coverage)) + 
#        
#        # primers
#        geom_segment(data = bed_primer, mapping = aes(x = chrom_start, xend = chrom_end, y = 0, yend = 0, color = as.factor(strand)), size = 8000) +
#        
#        geom_line() + 
#        #geom_area(alpha = 0.3) + # time intensive
#        
#        
#        facet_grid(sample~.) + 
#        theme_minimal() +
#        theme(legend.position = "none") +
#        xlim(0, coverage$position %>% max) + 
#        labs(); #p
#    
#    lst[[sample_name]] = plot_part
#    
#}
#
#lst[["scale"]] = b
#
#
#al = align_plots(
#    plotlist = lst,
#    align = "v",
#    axis = "lr",
#    greedy = TRUE
#) 
#
#plot_grid(plotlist = al, ncol = 1)
#
#
#
#ggsave("plots/cov.png", height = 10, width = 30)    
#
#
#