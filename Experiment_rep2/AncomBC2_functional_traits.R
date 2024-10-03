# ANCOM-BC2 functional traits
## CC

```{r phyloseq_ele, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_CC <- sample_metadata%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Control") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2") %>% 
  sample_data()

count_phy_CC <- GIFTs_elements_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(samples_CC)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

TAX <- uniqueGIFT_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_element")%>%
  as.matrix()%>%
  tax_table()

physeq_function_CC = phyloseq(count_phy_CC, TAX, samples_CC)
```

```{r phyloseq_ele_CC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_element_CC = ancombc2(data = physeq_function_CC, 
                                        assay_name = "counts", 
                                        tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                        fix_formula = "time_point", #fixed variable(s)
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
tax <- data.frame(physeq_function_CC@tax_table) %>%
  rownames_to_column(., "taxon")
```

```{r ancom_rand_res_elem1, echo=FALSE, comment="", message=FALSE, warning=FALSE}
colnames(ancom_rand_output_element_CC$res) <- gsub("-", "_", colnames(ancom_rand_output_element_CC$res))

ancombc_rand_table_CC <- ancom_rand_output_element_CC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_elem1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE}
ancombc_rand_table_CC%>%
  #      mutate(Function=factor(Function,levels=ancombc_rand_table$Function)) %>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
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

#### Functional level
```{r phylo_func1, comment="", echo=FALSE, message=FALSE, warning=FALSE}
sample_info_tab_phyCC <- accli_post2%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Control") %>%
  sample_data()

count_phyCC <- GIFTs_functions_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(sample_info_tab_phyCC)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

unique_funct_db<- GIFT_db[c(3,4,5)] %>% 
  distinct(Code_function, .keep_all = TRUE)

TAX <- unique_funct_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_function")%>%
  as.matrix()%>%
  tax_table()

physeq_functional_filtered_CC = phyloseq(count_phyCC, TAX, sample_info_tab_phyCC)
physeq_functional_filtered_clr_CC <- microbiome::transform(physeq_functional_filtered_CC, 'clr')  
```


```{r ancom_rand_func_prewild, comment="", echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_function_CC = ancombc2(data = physeq_functional_filtered_CC, 
                                         assay_name = "counts", 
                                         tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                         fix_formula = "time_point", #fixed variable(s)
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

```{r table_func1, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_functional_filtered_clr_CC@tax_table) %>%
  rownames_to_column(., "taxon")

colnames(ancom_rand_output_function_CC$res) <- gsub("-", "_", colnames(ancom_rand_output_function_CC$res))

ancombc_rand_func_CC <- ancom_rand_output_function_CC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon") %>%
  dplyr::arrange(lfc_time_point6_Post_FMT2)
```

```{r ancom_rand_res_funct1, comment="", echo=FALSE, message=FALSE, warning=FALSE}
ancombc_rand_table_func_CC <- ancom_rand_output_function_CC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_funct, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=4, fig.width=8, fig.fullwidth=TRUE}
ancombc_rand_table_func_CC%>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
  geom_col() +
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

## CI

```{r phyloseq_ele_1, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_CI <- sample_metadata%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Treatment") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2") %>% 
  sample_data()

count_phy_CI <- GIFTs_elements_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(samples_CI)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

TAX <- uniqueGIFT_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_element")%>%
  as.matrix()%>%
  tax_table()

physeq_function_CI = phyloseq(count_phy_CI, TAX, samples_CI)
```

```{r phyloseq_ele_CI, comment="", echo=FALSE, message=FALSE, warning=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_element_CI = ancombc2(data = physeq_function_CI, 
                                        assay_name = "counts", 
                                        tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                        fix_formula = "time_point", #fixed variable(s)
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

```{r table_ele1CI, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_function_CI@tax_table) %>%
  rownames_to_column(., "taxon")
```

```{r ancom_rand_res_elem1CI, echo=FALSE, comment="", message=FALSE, warning=FALSE}
colnames(ancom_rand_output_element_CI$res) <- gsub("-", "_", colnames(ancom_rand_output_element_CI$res))

ancombc_rand_table_CI <- ancom_rand_output_element_CI$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_elem1CI, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE}
ancombc_rand_table_CI%>%
  #      mutate(Function=factor(Function,levels=ancombc_rand_table$Function)) %>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
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

#### Functional level
```{r phylo_func1CI, comment="", echo=FALSE, message=FALSE, warning=FALSE}
sample_info_tab_phyCI <- accli_post2%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Treatment") %>%
  sample_data()

count_phyCI<- GIFTs_functions_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(sample_info_tab_phyCI)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

unique_funct_db<- GIFT_db[c(3,4,5)] %>% 
  distinct(Code_function, .keep_all = TRUE)

TAX <- unique_funct_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_function")%>%
  as.matrix()%>%
  tax_table()

physeq_functional_filtered_CI = phyloseq(count_phyCI, TAX, sample_info_tab_phyCI)
physeq_functional_filtered_clr_CI <- microbiome::transform(physeq_functional_filtered_CI, 'clr')  
```


```{r ancom_rand_func_prewildCI, comment="", echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_function_CI = ancombc2(data = physeq_functional_filtered_CI, 
                                         assay_name = "counts", 
                                         tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                         fix_formula = "time_point", #fixed variable(s)
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

```{r table_func1CI, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_functional_filtered_clr_CI@tax_table) %>%
  rownames_to_column(., "taxon")

colnames(ancom_rand_output_function_CI$res) <- gsub("-", "_", colnames(ancom_rand_output_function_CI$res))

ancombc_rand_func_CI <- ancom_rand_output_function_CI$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon") %>%
  dplyr::arrange(lfc_time_point6_Post_FMT2)
```

```{r ancom_rand_res_funct1CI, comment="", echo=FALSE, message=FALSE, warning=FALSE}
ancombc_rand_table_func_CI <- ancom_rand_output_function_CI$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_functCI, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=4, fig.width=8, fig.fullwidth=TRUE}
ancombc_rand_table_func_CI%>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
  geom_col() +
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

## WC

```{r phyloseq_ele_1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_WC <- sample_metadata%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Hot_control") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2") %>% 
  sample_data()

count_phy_WC <- GIFTs_elements_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(samples_WC)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

TAX <- uniqueGIFT_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_element")%>%
  as.matrix()%>%
  tax_table()

physeq_function_WC = phyloseq(count_phy_WC, TAX, samples_WC)
```

```{r phyloseq_ele_WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_element_WC = ancombc2(data = physeq_function_WC, 
                                        assay_name = "counts", 
                                        tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                        fix_formula = "time_point", #fixed variable(s)
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

```{r table_ele1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_function_WC@tax_table) %>%
  rownames_to_column(., "taxon")
```

```{r ancom_rand_res_elem1WC, echo=FALSE, comment="", message=FALSE, warning=FALSE}
colnames(ancom_rand_output_element_WC$res) <- gsub("-", "_", colnames(ancom_rand_output_element_WC$res))

ancombc_rand_table_WC <- ancom_rand_output_element_WC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_elem1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE}
ancombc_rand_table_WC%>%
  #      mutate(Function=factor(Function,levels=ancombc_rand_table$Function)) %>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
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

#### Functional level
```{r phylo_func1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
sample_info_tab_phyWC <- accli_post2%>%
  column_to_rownames(var="Tube_code")%>%
  filter(type=="Hot_control") %>%
  sample_data()

count_phyWC<- GIFTs_functions_community %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(sample_info_tab_phyWC)) %>%
  column_to_rownames("sample") %>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% 
  t() %>%
  otu_table(., taxa_are_rows=T)

unique_funct_db<- GIFT_db[c(3,4,5)] %>% 
  distinct(Code_function, .keep_all = TRUE)

TAX <- unique_funct_db%>%
  remove_rownames()%>%
  column_to_rownames(var="Code_function")%>%
  as.matrix()%>%
  tax_table()

physeq_functional_filtered_WC = phyloseq(count_phyWC, TAX, sample_info_tab_phyWC)
physeq_functional_filtered_clr_WC <- microbiome::transform(physeq_functional_filtered_WC, 'clr')  
```


```{r ancom_rand_func_prewildWC, comment="", echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
set.seed(1234) #set seed for reproducibility
ancom_rand_output_function_WC = ancombc2(data = physeq_functional_filtered_WC, 
                                         assay_name = "counts", 
                                         tax_level = NULL, #change to agglomerate analysis to a higher taxonomic range
                                         fix_formula = "time_point", #fixed variable(s)
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

```{r table_func1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
tax <- data.frame(physeq_functional_filtered_clr_WC@tax_table) %>%
  rownames_to_column(., "taxon")

colnames(ancom_rand_output_function_WC$res) <- gsub("-", "_", colnames(ancom_rand_output_function_WC$res))

ancombc_rand_func_WC <- ancom_rand_output_function_WC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon") %>%
  dplyr::arrange(lfc_time_point6_Post_FMT2)
```

```{r ancom_rand_res_funct1WC, comment="", echo=FALSE, message=FALSE, warning=FALSE}
ancombc_rand_table_func_WC <- ancom_rand_output_function_WC$res %>%
  dplyr::select(taxon, lfc_time_point6_Post_FMT2, p_time_point6_Post_FMT2) %>%
  filter(p_time_point6_Post_FMT2 < 0.05) %>%
  dplyr::arrange(p_time_point6_Post_FMT2) %>%
  merge(., tax, by="taxon")
```

```{r ancombc_rand_plot_functWC, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=4, fig.width=8, fig.fullwidth=TRUE}
ancombc_rand_table_func_WC%>%
  mutate(Color = ifelse(lfc_time_point6_Post_FMT2 <0, "Acclimation","Post2")) %>%
  ggplot(aes(x=forcats::fct_reorder(Function,lfc_time_point6_Post_FMT2), y=lfc_time_point6_Post_FMT2, fill=Color)) + 
  geom_col() +
  scale_fill_manual(values=c("#6b7398","#e5bd5b")) + 
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