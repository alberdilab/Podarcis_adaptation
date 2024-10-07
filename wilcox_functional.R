## Wilcoxon comparison

### Community elements differences: in CC accli vs post2 ####
```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_CC <- sample_metadata%>%
  filter(type=="Control") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2")%>% 
  mutate(time_point = recode(time_point,
                        "1_Acclimation" = "Acclimation", 
                        "6_Post-FMT2" = "Post2"))

element_gift_CC <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(samples_CC[c(1,10)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_CC <- element_gift_CC %>%
  pivot_longer(-c(Tube_code,time_point), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ time_point)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_value < 0.05)%>% #note that p_value is used instead of p_adjust
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_CC <- element_gift_CC  %>% 
  dplyr::select(-c(time_point))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_CC <- subset(element_gift_t_CC, trait %in% significant_elements_CC$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CC[c(1,10)], by = join_by(Tube_code == Tube_code))

difference_table_CC <- element_gift_filt_CC %>%
  select(-Tube_code) %>%
  group_by(time_point) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Acclimation-Post2)%>% 
  mutate(group_color = ifelse(Difference <0, "Post2","Acclimation")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_CC <- element_gift_filt_CC %>% 
  select(-time_point) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(samples_CC, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(time_point, elements) %>%
  summarise(mean=mean(abundance))

log_fold_CC <- means_gift_CC %>%
  group_by(elements) %>%
  summarise(
    logfc_accli_post2 = log2(mean[time_point == "Acclimation"] / mean[time_point == "Post2"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CC, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_value), color=group_color)) + #note that p_value was used instead of p_adjust
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_CC <- difference_table_CC %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_CC<-unique(code_function_CC$Code_function)

gift_colors_CC <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_CC)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_CC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CC, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_CC, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_accli_post2, y = -log(p_value), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_CC$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_CC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CC, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_CC, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_accli_post2, y = -log(p_value), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_CC$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CC, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_accli_post2, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CC %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```
```{r}
difference_table_CC %>% 
  filter(group_color=="Acclimation")
```
```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
sugar_function <- GIFT_db %>% 
  filter(Code_function=="D03") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(sugar_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
poly_function <- GIFT_db %>% 
  filter(Code_function=="D02") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(poly_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r elements_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_names_CC <- element_gift_filt_CC%>%
  select(-time_point)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CC[c(1,10)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_CC)[2:16] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_CC, aes(x=time_point, y=.data[[i]], color = time_point, fill=time_point)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}

```

### Community elements differences: in CI accli vs post2 ####

```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_CI <- sample_metadata%>%
  filter(type=="Treatment") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2")%>% 
  mutate(time_point = recode(time_point,
                             "1_Acclimation" = "Acclimation", 
                             "6_Post-FMT2" = "Post2"))

element_gift_CI <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(samples_CI[c(1,10)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_CI <- element_gift_CI %>%
  pivot_longer(-c(Tube_code,time_point), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ time_point)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>%
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_CI <- element_gift_CI  %>% 
  dplyr::select(-c(time_point))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_CI <- subset(element_gift_t_CI, trait %in% significant_elements_CI$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CI[c(1,10)], by = join_by(Tube_code == Tube_code))

difference_table_CI <- element_gift_filt_CI %>%
  select(-Tube_code) %>%
  group_by(time_point) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Acclimation-Post2)%>% 
  mutate(group_color = ifelse(Difference <0, "Post2","Acclimation")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_CI <- element_gift_filt_CI %>% 
  select(-time_point) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(samples_CI, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(time_point, elements) %>%
  summarise(mean=mean(abundance))

log_fold_CI <- means_gift_CI %>%
  group_by(elements) %>%
  summarise(
    logfc_acclimation_post2 = log2(mean[time_point == "Acclimation"] / mean[time_point == "Post2"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CI %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CI, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_adjust), color=group_color)) + 
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_CI <- difference_table_CI %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_CI<-unique(code_function_CI$Code_function)

gift_colors_CI <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_CI)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_CI %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CI, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CI, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_CI, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_acclimation_post2, y = -log(p_adjust), color=legend, size=abs(Difference))) +
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_CI$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_CI %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CI, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CI, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_CI, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_acclimation_post2, y = -log(p_adjust), color=legend, size=abs(Difference))) + 
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_CI$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CI %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_CI, by=join_by(Elements==trait)) %>%
  left_join(log_fold_CI, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_acclimation_post2, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=FALSE}
difference_table_CI %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```
```{r}
difference_table_CI %>% 
  filter(group_color=="Acclimation")
```
```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
sugar_function <- GIFT_db %>% 
  filter(Code_function=="D03") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(sugar_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
poly_function <- GIFT_db %>% 
  filter(Code_function=="D02") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(poly_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r elements_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_names_CI <- element_gift_filt_CI%>%
  select(-time_point)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CI[c(1,10)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_CI)[2:18] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_CI, aes(x=time_point, y=.data[[i]], color = time_point, fill=time_point)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}

```

### Community elements differences: in WC accli vs post2 ####
```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_WC <- sample_metadata%>%
  filter(type=="Hot_control") %>% 
  filter(time_point == "1_Acclimation"|time_point == "6_Post-FMT2")%>% 
  mutate(time_point = recode(time_point,
                             "1_Acclimation" = "Acclimation", 
                             "6_Post-FMT2" = "Post2"))

element_gift_WC <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(samples_WC[c(1,10)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_WC <- element_gift_WC %>%
  pivot_longer(-c(Tube_code,time_point), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ time_point)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_value < 0.05)%>% #note that p_value is used instead of p_adjust
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_WC <- element_gift_WC  %>% 
  dplyr::select(-c(time_point))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_WC <- subset(element_gift_t_WC, trait %in% significant_elements_WC$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_WC[c(1,10)], by = join_by(Tube_code == Tube_code))

difference_table_WC <- element_gift_filt_WC %>%
  select(-Tube_code) %>%
  group_by(time_point) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Acclimation-Post2)%>% 
  mutate(group_color = ifelse(Difference <0, "Post2","Acclimation")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_WC <- element_gift_filt_WC %>% 
  select(-time_point) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(samples_WC, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(time_point, elements) %>%
  summarise(mean=mean(abundance))

log_fold_WC <- means_gift_WC %>%
  group_by(elements) %>%
  summarise(
    logfc_acclimation_post2 = log2(mean[time_point == "Acclimation"] / mean[time_point == "Post2"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_WC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_WC, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_value), color=group_color)) + #note that p_value was used instead of p_adjust
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_WC <- difference_table_WC %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_WC<-unique(code_function_WC$Code_function)

gift_colors_WC <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_WC)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_WC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_WC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_WC, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_WC, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_acclimation_post2, y = -log(p_value), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_WC$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_WC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_WC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_WC, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_WC, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_acclimation_post2, y = -log(p_value), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_WC$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_WC %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_WC, by=join_by(Elements==trait)) %>%
  left_join(log_fold_WC, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_acclimation_post2, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=FALSE}
difference_table_WC %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```
```{r}
difference_table_WC %>% 
  filter(group_color=="Acclimation")
```
```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
sugar_function <- GIFT_db %>% 
  filter(Code_function=="D03") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(sugar_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
poly_function <- GIFT_db %>% 
  filter(Code_function=="D02") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(poly_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r elements_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_names_WC <- element_gift_filt_WC%>%
  select(-time_point)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_WC[c(1,10)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_WC)[2:5] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_WC, aes(x=time_point, y=.data[[i]], color = time_point, fill=time_point)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}

```

### Comparison of both population in wild samples ####

```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_wild<- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(sample_metadata_wild[c(1,4)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_wild <- element_gift_accli %>%
  pivot_longer(-c(Tube_code,Population), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ Population)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>%
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_wild <- element_gift_wild  %>% 
  dplyr::select(-c(Population))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_wild <- subset(element_gift_t_wild, trait %in% significant_elements_wild$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_wild[c(1,4)], by = join_by(Tube_code == Tube_code))

difference_table_wild <- element_gift_filt_wild %>%
  select(-Tube_code) %>%
  group_by(Population) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Cold_wet-Hot_dry)%>% 
  mutate(group_color = ifelse(Difference <0, "Hot","Cold")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_wild <- element_gift_filt_wild %>% 
  select(-Population) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(sample_metadata_wild, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(Population, elements) %>%
  summarise(mean=mean(abundance))

log_fold_wild <- means_gift_wild %>%
  group_by(elements) %>%
  summarise(
    logfc_Cold_hot = log2(mean[Population == "Cold_wet"] / mean[Population == "Hot_dry"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_wild %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_wild, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_adjust), color=group_color)) + 
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_wild <- difference_table_wild %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_wild<-unique(code_function_wild$Code_function)

gift_colors_wild <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_wild)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_wild %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_wild, by=join_by(Elements==trait)) %>%
  left_join(log_fold_wild, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_wild, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_Cold_hot, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_wild$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_wild %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_wild, by=join_by(Elements==trait)) %>%
  left_join(log_fold_wild, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_wild, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_Cold_hot, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_wild$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_wild %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_wild, by=join_by(Elements==trait)) %>%
  left_join(log_fold_wild, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_Cold_hot, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=FALSE}
difference_table_wild %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```
```{r}
difference_table_WC %>% 
  filter(group_color=="Acclimation")
```
```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
sugar_function <- GIFT_db %>% 
  filter(Code_function=="D03") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(sugar_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
poly_function <- GIFT_db %>% 
  filter(Code_function=="D02") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(poly_function)) %>% 
  rownames_to_column(., "sample") %>% 
  left_join(sample_metadata[c(1,21)], by=join_by(sample==sample)) %>% 
  select(-sample) %>% 
  pivot_longer(!diet,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=diet, y=gift, color=diet)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.5,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=diet_colors)+
  scale_fill_manual(values=diet_colors) +
  scale_x_discrete(labels=c("Wild" = "Wild", "Pre_grass" = "Captive-born")) +
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r elements_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_names_wild <- element_gift_filt_wild%>%
  select(-Population)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_wild[c(1,4)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_wild)[2:18] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_wild, aes(x=Population, y=.data[[i]], color = Population, fill=Population)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}

```

#Butiryc acid biosynthesis

lipid_function <- GIFT_db %>% 
  filter(Code_function=="D01") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(lipid_function)) %>% 
  rownames_to_column(., "Tube_code") %>% 
  left_join(samples_WC[c(1,10)], by=join_by(Tube_code==Tube_code)) %>%
  filter(time_point!="NA") %>% 
  select(-Tube_code) %>% 
  pivot_longer(!time_point,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=time_point, y=gift, color=time_point, fill=time_point)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.3,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  scale_fill_manual(values=c("#e5bd5b50", "#6b739850"))+
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )


### Comparison of both population in acclimation samples ####

```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_accli<- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(sample_metadata_accli[c(1,4)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_accli <- element_gift_accli %>%
  pivot_longer(-c(Tube_code,Population), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ Population)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>%
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_accli <- element_gift_accli  %>% 
  dplyr::select(-c(Population))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_accli <- subset(element_gift_t_accli, trait %in% significant_elements_accli$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_accli[c(1,4)], by = join_by(Tube_code == Tube_code))

difference_table_accli <- element_gift_filt_accli %>%
  select(-Tube_code) %>%
  group_by(Population) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Cold_wet-Hot_dry)%>% 
  mutate(group_color = ifelse(Difference <0, "Hot","Cold")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_accli <- element_gift_filt_accli %>% 
  select(-Population) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(sample_metadata_accli, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(Population, elements) %>%
  summarise(mean=mean(abundance))

log_fold_accli <- means_gift_accli %>%
  group_by(elements) %>%
  summarise(
    logfc_Cold_hot = log2(mean[Population == "Cold_wet"] / mean[Population == "Hot_dry"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_accli %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_accli, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_adjust), color=group_color)) + 
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_accli <- difference_table_accli %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_accli<-unique(code_function_accli$Code_function)

gift_colors_accli <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_accli)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_accli %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_accli, by=join_by(Elements==trait)) %>%
  left_join(log_fold_accli, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_accli, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_Cold_hot, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_accli$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_accli %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_accli, by=join_by(Elements==trait)) %>%
  left_join(log_fold_accli, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_accli, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_Cold_hot, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_accli$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_accli %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_accli, by=join_by(Elements==trait)) %>%
  left_join(log_fold_accli, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_Cold_hot, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=FALSE}
difference_table_accli %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```

element_gift_names_CC <- element_gift_filt_CC%>%
  select(-time_point)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CC[c(1,10)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_CC)[2:16] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_CC, aes(x=time_point, y=.data[[i]], color = time_point, fill=time_point)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}


### Comparison of post2 samples of CC vs CI ####

```{r comunity_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
samples_tm6 <- sample_metadata%>%
  filter(type!="Hot_control") %>% 
  filter(time_point == "6_Post-FMT2")%>% 
  mutate(time_point = recode(time_point,
                             "6_Post-FMT2" = "Post2"))

element_gift_tm6 <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(samples_tm6[c(1,7)], by="Tube_code") 
```

```{r commun_wilcox_elem, comment="", echo=FALSE, message=FALSE, warning=FALSE}
uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)

significant_elements_tm6 <- element_gift_tm6 %>%
  pivot_longer(-c(Tube_code,type), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ type)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>% #note that p_value is used instead of p_adjust
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)
# %>%
#   left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t_tm6 <- element_gift_tm6  %>% 
  dplyr::select(-c(type))  %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt_tm6 <- subset(element_gift_t_tm6, trait %in% significant_elements_tm6$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_tm6[c(1,7)], by = join_by(Tube_code == Tube_code))

difference_table_tm6 <- element_gift_filt_tm6 %>%
  select(-Tube_code) %>%
  group_by(type) %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(., "Elements") %>%
  left_join(.,uniqueGIFT_db[c(1,3,4)],by = join_by(Elements == Code_element)) %>% 
  arrange(Function) %>% 
  mutate(Difference=Control-Treatment)%>% 
  mutate(group_color = ifelse(Difference <0, "Treatment","Control")) 
```

```{r log_fold_calc, comment="", echo=FALSE, message=FALSE, warning=FALSE}
means_gift_tm6 <- element_gift_filt_tm6 %>% 
  select(-type) %>% 
  pivot_longer(!Tube_code, names_to = "elements", values_to = "abundance") %>% 
  left_join(samples_tm6, by=join_by(Tube_code==Tube_code)) %>% 
  group_by(type, elements) %>%
  summarise(mean=mean(abundance))

log_fold_tm6 <- means_gift_tm6 %>%
  group_by(elements) %>%
  summarise(
    logfc_control_treatment = log2(mean[type == "Control"] / mean[type == "Treatment"])
  )
```


```{r plot1, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_tm6 %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_tm6, by=join_by(Elements==trait)) %>% 
  ggplot(., aes(x = Difference, y = -log(p_adjust), color=group_color)) + #note that p_value was used instead of p_adjust
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2)
```

```{r plot3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=5, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
uniqueGIFT <- unique(GIFT_db[c(2,3,4,5,6)])

code_function_tm6 <- difference_table_tm6 %>%
  left_join(uniqueGIFT[c(1:3)], by=join_by(Elements==Code_element))

unique_codes_tm6<-unique(code_function_tm6$Code_function)

gift_colors_tm6 <- read_tsv("data/gift_colors.tsv") %>% 
  filter(Code_function %in% unique_codes_tm6)%>% 
  mutate(legend=str_c(Code_function," - ",Function))

code_function_tm6 %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_tm6, by=join_by(Elements==trait)) %>%
  left_join(log_fold_tm6, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_tm6, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_control_treatment, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_point()+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_tm6$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value")

# +
#   geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)

```
```{r final_plot_3, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=12, fig.fullwidth=TRUE}
code_function_tm6 %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_tm6, by=join_by(Elements==trait)) %>%
  left_join(log_fold_tm6, by=join_by(Elements==elements)) %>% 
  left_join(gift_colors_tm6, by=join_by(Code_function==Code_function)) %>% 
  ggplot(., aes(x = logfc_control_treatment, y = -log(p_adjust), color=legend, size=abs(Difference))) + #note that p_value was used instead of p_adjust
  geom_jitter(width = 0.2, height = 0.2)+
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors_tm6$Color)+
  #xlim(c(-10,4)) +
  theme_classic()+
  labs(size="Mean difference (abs)", color="Functional trait")+
  labs(x = "Log-fold change", y="-Log adjusted p-value") +
  geom_text_repel(aes(label = Element), min.segment.length = 0.4, size=2.5, max.overlaps = Inf)
```

```{r plot4, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=8, fig.width=10, fig.fullwidth=TRUE, eval=FALSE}
difference_table_tm6 %>%
  #  mutate(Difference_abs = abs(Difference)) %>% 
  left_join(significant_elements_tm6, by=join_by(Elements==trait)) %>%
  left_join(log_fold_tm6, by=join_by(Elements==elements)) %>% 
  ggplot(., aes(x = Difference, y = logfc_control_treatment, color=group_color)) +
  geom_point(size=2, show.legend = FALSE) + 
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  #xlim(c(-10,4)) +
  theme_classic()+
  geom_text_repel(aes(label = Elements), size=2, max.overlaps = Inf)
```

```{r commun_wilcox_elem_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE, eval=TRUE}
difference_table_tm6 %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
  #  geom_point(size=4) + 
  scale_fill_manual(values=c("#e5bd5b", "#6b7398")) + 
  geom_hline(yintercept=0) + 
  coord_flip()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "right", 
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "grey"))+
  xlab("Function") + 
  ylab("Mean difference")
```
```{r}
difference_table_CC %>% 
  filter(group_color=="Acclimation")
```
```{r sugar_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
sugar_function <- GIFT_db %>% 
  filter(Code_function=="D03") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(sugar_function)) %>% 
  rownames_to_column(., "Tube_code") %>% 
  left_join(samples_CC[c(1,10)], by=join_by(Tube_code==Tube_code)) %>%
  filter(time_point!="NA") %>% 
  select(-Tube_code) %>% 
  pivot_longer(!time_point,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=time_point, y=gift, color=time_point, fill=time_point)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.3,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  scale_fill_manual(values=c("#e5bd5b50", "#6b739850"))+
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r poly_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
poly_function <- GIFT_db %>% 
  filter(Code_function=="D02") %>% 
  select(Code_element) %>% 
  pull()

GIFTs_elements_community %>% 
  as.data.frame() %>% 
  select(any_of(poly_function)) %>% 
  rownames_to_column(., "Tube_code") %>% 
  left_join(samples_CC[c(1,10)], by=join_by(Tube_code==Tube_code)) %>%
  filter(time_point!="NA") %>% 
  select(-Tube_code) %>% 
  pivot_longer(!time_point,names_to="trait",values_to="gift") %>% 
  left_join(uniqueGIFT_db, by=join_by(trait==Code_element)) %>% 
  ggplot(aes(x=time_point, y=gift, color=time_point, fill=time_point)) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  geom_boxplot(width = 0.5, alpha=0.3,outlier.shape = NA, show.legend = FALSE) +
  scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
  scale_fill_manual(values=c("#e5bd5b50", "#6b739850"))+
  facet_grid(.~Element,  scales="free_x")+
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text = element_text(size=10, angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  )
```

```{r elements_plot, comment="", echo=FALSE, message=FALSE, warning=FALSE}
element_gift_names_CC <- element_gift_filt_CC%>%
  select(-time_point)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  select(-Elements)%>%
  select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., samples_CC[c(1,10)], by = join_by(Tube_code == Tube_code))


colNames <- names(element_gift_names_CC)[2:16] #check the column that has the last code
for(i in colNames){
  plt <- ggplot(element_gift_names_CC, aes(x=time_point, y=.data[[i]], color = time_point, fill=time_point)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    scale_color_manual(values=c("#e5bd5b", "#6b7398"))+
    scale_fill_manual(values=c("#e5bd5b", "#6b7398"))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}

```
