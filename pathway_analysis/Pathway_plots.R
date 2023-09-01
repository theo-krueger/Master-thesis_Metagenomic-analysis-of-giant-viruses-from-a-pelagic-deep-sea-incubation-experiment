### setup
library(optparse)

parser = OptionParser()

parser <- add_option(parser,
                     c("-i", "--input-file"), 
                     default = NULL,
                     help = "Metabolic pathway summary file")
parser <- add_option(parser,
                     c("-o", "--output-dir"), 
                     type="character", 
                     default=NULL,
                     help="Output directory [default = 'input_directory'/figures]")

args <- parse_args(parser)

#print(args)

if (is.null(args$`input-file`)){
  print_help(parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

suppressPackageStartupMessages(library(tidyverse))

#logical arguments
`%notin%` <- Negate(`%in%`)

date = format(Sys.Date(), "%Y%m%d")

# paths 
path_input <-args$`input-file`[1]

if (is.null(args$`output-dir`)){
  path_output <- dirname(path_input)
  path_output_figures <- paste0(path_output, "/figures/")
  
  if (!dir.exists(path_output_figures)){
    dir.create(path_output_figures)
  } 
} else {
  path_output_figures <- args$`output-dir`[1]
}

# mode
save <- TRUE

# load
df_pathways <- read.table(path_input, sep = "\t", header = TRUE) %>%
  complete(genome_name, module_name)

module_lookup <- df_pathways %>%
  select(module_name, module_category, module_subcategory)%>%
  distinct() %>%
  na.omit()

df_pathways_filled <- df_pathways %>%
  select(-module_category, -module_subcategory) %>%
  left_join(module_lookup, by="module_name")
  # mutate(module_completeness = case_when(is.na(module_completeness) ~ 0,
  #                                        TRUE ~ module_completeness))

plot_list = unique(df_pathways_filled$module_category)

for (plot_type in plot_list){
  
  module_name <- gsub(pattern = " ", replacement = "_", x = plot_type)
  
  print(paste0("Processing: ", plot_type))
  
  df_sub <- df_pathways_filled %>%
    filter(module_category == plot_type) 
  
  # plot
plot <- ggplot(data = df_sub)+
  geom_point(aes(x = genome_name, 
                 y = module_name, 
                 alpha = module_completeness,
                 group = module_subcategory),
             size = 5) + 
  # geom_hline(yintercept = seq(1.5, length(unique(df_sub$module_name))-0.5, 1), 
  #            lwd = 0.2)+
  # geom_vline(xintercept = seq(1.5, length(unique(df_sub$genome_name))-0.5, 1), 
  #            lwd = 0.2)+
  facet_grid(module_subcategory~.,
             switch = "y",
             scales = "free",
             space = "free",
             labeller = label_wrap_gen())+
  # annotate("segment",x=-Inf,xend=-Inf,y=-Inf,yend=Inf,color="black",lwd=1)+
  scale_y_discrete(labels = label_wrap_gen())+
  scale_alpha(limits = c(0,1), 
              breaks = seq(0,1,0.2), 
              na.value = 0)+
  labs(x = "Bin", 
       y = "KEGG Module Name", 
       alpha = "Completeness", 
       title = plot_type)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        strip.text = element_text(margin = margin(2, 0.25, 2, 0.25, "cm")),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  
  if (save) {
    ggsave(plot = plot, 
           filename = paste0(path_output_figures, date, "_", "pathway_", module_name, ".pdf"), 
           device = "pdf",
           width = 10,
           dpi = 300)
  }
}

