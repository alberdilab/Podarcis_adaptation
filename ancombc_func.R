sample_info_tab_phy <- sample_metadata%>%
  filter(type=="Control") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2")%>% 
  mutate(time_point = recode(time_point,
                             "1_Acclimation" = "Acclimation", 
                             "6_Post-FMT2" = "Post2"))%>% 
  column_to_rownames("Tube_code")%>% 
  sample_data()

count_phy <- GIFTs_elements_community %>% 
  as.data.frame() %>%
  rownames_to_column("Tube_code") %>%
  filter(Tube_code %in% rownames(sample_info_tab_phy)) %>%
  column_to_rownames("Tube_code") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

uniqueGIFT_db<-uniqueGIFT_db[,-1]

TAX <- uniqueGIFT_db%>%
  remove_rownames%>%
  
  uniqueGIFT_db<-column_to_rownames(uniqueGIFT_db,"Code_element")%>%
  as.matrix()%>%
  tax_table()

physeq_function = phyloseq(count_phy, TAX, sample_info_tab_phy)

```

```{r ancom_rand_elem_capwild, comment="", echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_element_captive_wild = ancombc2(data = physeq_function, 
                                                  assay_name = "counts", 
                                                  tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                                  fix_formula = "region", #fixed variable(s)
                                                  p_adj_method = "holm", 
                                                  pseudo_sens = TRUE,
                                                  prv_cut = 0, 
                                                  lib_cut = 0, 
                                                  s0_perc = 0.05,
                                                  group = NULL, 
                                                  struc_zero = FALSE, 
                                                  neg_lb = FALSE,
                                                  alpha = 0.05, 
                                                  n_cl = 2, 
                                                  verbose = TRUE,
                                                  global = FALSE, 
                                                  pairwise = FALSE, 
                                                  dunnet = FALSE, 
                                                  trend = FALSE,
                                                  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
                                                  em_control = list(tol = 1e-5, max_iter = 100),
                                                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                                                  trend_control = NULL)

```

```{r table_ele1, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_function@tax_table) %>%
  rownames_to_column(., "taxon")
```

```{r ancom_rand_res_elem1, echo=FALSE, comment="", message=FALSE, warning=FALSE}
ancombc_rand_table <- ancom_rand_output_element_captive_wild$res %>%
  dplyr::select(taxon, lfc_regionNafarroa, p_regionNafarroa) %>%
  filter(p_regionNafarroa < 0.05) %>%
  dplyr::arrange(p_regionNafarroa) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_elem1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE}
ancombc_rand_table%>%
  #      mutate(Function=factor(Function,levels=ancombc_rand_table$Function)) %>%
  mutate(Color = ifelse(lfc_regionNafarroa <0, "Wild","Captive")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_regionNafarroa), y=lfc_regionNafarroa, fill=Color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10, face="bold.italic"),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Traits") + 
  ylab("log2FoldChange")
```