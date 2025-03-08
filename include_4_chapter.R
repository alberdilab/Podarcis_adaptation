# code to include in 4th chapter

### Do each population retain the wild GM after acclimation?
#### Cold
##### Alpha diversity
```{r alpha_div_plot_cold, comment="", message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(Population == "Cold_wet" & time_point %in% c("Acclimation", "Wild")) %>%
  mutate(time_point = factor(time_point, levels = c("Wild", "Acclimation"))) %>%
  mutate(metric = factor(metric, levels = c("richness", "neutral", "phylogenetic"))) %>%
  ggplot(aes(
    y = value,
    x = time_point,
    group = time_point,
    color = time_point,
    fill = time_point
  )) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(
    name = "Population",
    breaks = c("Wild", "Acclimation"),
    labels = c("Wild", "Acclimation"),
    values = c('#008080', "#d57d2c")
  ) +
  scale_fill_manual(
    name = "Population",
    breaks = c("Wild", "Acclimation"),
    labels = c("Wild", "Acclimation"),
    values = c('#00808050', "#d57d2c50")
  ) +
  facet_wrap(. ~ metric, scales = "free") +
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  ) +
  ylab("Alpha diversity")
```

```{r, wild_accli}
alpha_div_meta <- alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(Population=="Cold_wet") %>%
  filter(time_point =="Acclimation"| time_point == "Wild") 
```

***Richness***
  ```{r alpha_div_data_wild_accli, comment="", message=FALSE, warning=FALSE}
Modelq0GLMMNB <- MASS::glm.nb(richness ~ time_point, data = alpha_div_meta,trace=TRUE)
summary(Modelq0GLMMNB)
emmeans::emmeans(Modelq0GLMMNB, pairwise ~ time_point)
```

***Neutral***
  ```{r alpha_div_data_1_wild_accli, comment="", message=FALSE, warning=FALSE}
Modelq1n <- lm(formula = neutral ~ time_point, data = alpha_div_meta) 
summary(Modelq1n)
```

***Phylogenetic***
  ```{r alpha_div_data_2_wild_accli, comment="", message=FALSE, warning=FALSE}
Model_phylo <- lm(formula = phylogenetic ~ time_point, data = alpha_div_meta) 
summary(Model_phylo)
```

##### Beta diversity

```{r beta_div_neutral_accli_wild_cold, message=FALSE, warning=FALSE, comments=FALSE}
samples_to_keep_accli <- sample_metadata %>%
  filter(Population=="Cold_wet" & time_point %in% c("Acclimation", "Wild")) %>% 
  select(Tube_code) %>% 
  pull()
subset_meta_accli <- sample_metadata %>%
  filter(Population=="Cold_wet" & time_point %in% c("Acclimation", "Wild"))
```

***Neutral***
  ```{r beta_neutral_permanova_1_accli_cold, warning=FALSE, comments="", message=FALSE}
neutral_accli <- as.matrix(beta_q1n$S)
neutral_accli <- as.dist(neutral_accli[rownames(neutral_accli) %in% samples_to_keep_accli,
                                       colnames(neutral_accli) %in% samples_to_keep_accli])
betadisper(neutral_accli, subset_meta_accli$time_point) %>% permutest(., pairwise = TRUE)
```

```{r beta_neutral_permanova_accli, warning=FALSE, comments="", message=FALSE, results=TRUE}
adonis2(neutral_accli ~ time_point,
        data = subset_meta_accli %>% arrange(match(Tube_code,labels(neutral_accli))),
        permutations = 999,
        strata = subset_meta_accli %>% arrange(match(Tube_code,labels(neutral_accli))) %>% pull(individual)) %>%
  broom::tidy() %>%
  tt()
```

##### Functional differences
***MCI***
  ```{r gitfs_elements_accli_cold, warning=FALSE, comments="", message=FALSE}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(Population == "Cold_wet" & time_point %in% c("Acclimation", "Wild")) %>%
  group_by(time_point) %>%
  summarise(MCI = mean(value), sd = sd(value))
```
```{r gitfs_functional_accli_treat, warning=FALSE, comments="", message=FALSE}
MCI_element_accli <- GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(type == "Treatment" & time_point %in% c("Acclimation", "Wild"))

shapiro.test(MCI_element_accli$value)
var.test(value ~ time_point, data = MCI_element_accli)
t.test(value ~ time_point, data=MCI_element_accli, var.equal = TRUE)
```

***Differences in functional pathways***
  ```{r comunity_elem_fun, comment="", message=FALSE, warning=FALSE}
element_gift_wild <- element_gift %>% 
  filter(type %in% c("Control", "Treatment") & time_point %in% c("Wild", "Acclimation")) %>% 
  select(-time_point,-type) %>%
  filter(rowSums(across(where(is.numeric), ~ . != 0)) > 0) %>%
  select(Tube_code, where(~ is.numeric(.) && sum(.) > 0)) %>%
  left_join(sample_metadata[c(1,10)], by = "Tube_code")
```
```{r commun_wilcox_elem_wild_acc, comment="", message=FALSE, warning=FALSE}
significant_elements_wild <- element_gift_wild %>%
  pivot_longer(-c(Tube_code,time_point), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ time_point)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>%
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)

element_gift_sig_wild <- element_gift_wild %>%
  select(Tube_code, all_of(intersect(
    significant_elements_wild$trait,
    colnames(element_gift_wild)
  ))) %>%
  left_join(., sample_metadata[c(1, 10)], by = join_by(Tube_code == Tube_code))

difference_table_wild <- element_gift_sig_wild %>%
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
  mutate(Difference=Wild-Acclimation)%>% 
  mutate(group_color = ifelse(Difference <0, "Acclimation","Wild")) 
difference_table_wild
```

#### Warm
##### Alpha diversity
```{r alpha_div_plot_warm, comment="", message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(type=="Hot_control" & time_point %in% c("Acclimation", "Wild")) %>%
  mutate(time_point = factor(time_point, levels = c("Wild", "Acclimation"))) %>%
  mutate(metric = factor(metric, levels = c("richness", "neutral", "phylogenetic"))) %>%
  ggplot(aes(
    y = value,
    x = time_point,
    group = time_point,
    color = time_point,
    fill = time_point
  )) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(
    name = "Population",
    breaks = c("Wild", "Acclimation"),
    labels = c("Wild", "Acclimation"),
    values = c('#008080', "#d57d2c")
  ) +
  scale_fill_manual(
    name = "Population",
    breaks = c("Wild", "Acclimation"),
    labels = c("Wild", "Acclimation"),
    values = c('#00808050', "#d57d2c50")
  ) +
  facet_wrap(. ~ metric, scales = "free") +
  coord_cartesian(xlim = c(1, NA)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  ) +
  ylab("Alpha diversity")
```

```{r, warm_wild_accli}
alpha_div_meta <- alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(Population=="Hot_dry") %>%
  filter(time_point =="Acclimation"| time_point == "Wild") 
```

***Richness***
  ```{r alpha_div_data_warm_wild_accli, comment="", message=FALSE, warning=FALSE}
Modelq0GLMMNB <- MASS::glm.nb(richness ~ time_point, data = alpha_div_meta,trace=TRUE)
summary(Modelq0GLMMNB)
emmeans::emmeans(Modelq0GLMMNB, pairwise ~ time_point)
```

***Neutral***
  ```{r alpha_div_data_1_warm_wild_accli, comment="", message=FALSE, warning=FALSE}
Modelq1n <- lm(formula = neutral ~ time_point, data = alpha_div_meta) 
summary(Modelq1n)
```

***Phylogenetic***
  ```{r alpha_div_data_2_warm_wild_accli, comment="", message=FALSE, warning=FALSE}
Model_phylo <- lm(formula = phylogenetic ~ time_point, data = alpha_div_meta) 
summary(Model_phylo)
```

##### Beta diversity
```{r beta_div_neutral_accli_wild_warm, message=FALSE, warning=FALSE, comments=FALSE}
samples_to_keep_accli <- sample_metadata %>%
  filter(Population=="Hot_dry" & time_point %in% c("Acclimation", "Wild")) %>% 
  select(Tube_code) %>% 
  pull()
subset_meta_accli <- sample_metadata %>%
  filter(Population=="Hot_dry" & time_point %in% c("Acclimation", "Wild"))
```

***Neutral***
  ```{r beta_neutral_permanova_1_accli, warning=FALSE, comments="", message=FALSE}
neutral_accli <- as.matrix(beta_q1n$S)
neutral_accli <- as.dist(neutral_accli[rownames(neutral_accli) %in% samples_to_keep_accli,
                                       colnames(neutral_accli) %in% samples_to_keep_accli])
betadisper(neutral_accli, subset_meta_accli$time_point) %>% permutest(., pairwise = TRUE)
```

```{r beta_neutral_permanova_accli_hot, warning=FALSE, comments="", message=FALSE, results=TRUE}
adonis2(neutral_accli ~ time_point,
        data = subset_meta_accli %>% arrange(match(Tube_code,labels(neutral_accli))),
        permutations = 999,
        strata = subset_meta_accli %>% arrange(match(Tube_code,labels(neutral_accli))) %>% pull(individual)) %>%
  broom::tidy() %>%
  tt()
```

##### Functional differences
***MCI***
  ```{r gitfs_elements_accli_hot, warning=FALSE, comments="", message=FALSE}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(type == "Hot_control" & time_point %in% c("Acclimation", "Wild")) %>%
  group_by(time_point) %>%
  summarise(MCI = mean(value), sd = sd(value))
```
```{r gitfs_functional_accli_hot, warning=FALSE, comments="", message=FALSE}
MCI_element_accli <- GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(type == "Hot_control" & time_point %in% c("Acclimation", "Wild"))


shapiro.test(MCI_element_accli$value)
var.test(value ~ time_point, data = MCI_element_accli)
t.test(value ~ time_point, data=MCI_element_accli, var.equal = TRUE)
```

***Differences in functional pathways***
  ```{r comunity_elem_warm, comment="", message=FALSE, warning=FALSE}
element_gift_warm <- element_gift %>% 
  filter(type=="Hot_control" & time_point %in% c("Wild", "Acclimation")) %>% 
  select(-time_point,-type) %>%
  filter(rowSums(across(where(is.numeric), ~ . != 0)) > 0) %>%
  select(Tube_code, where(~ is.numeric(.) && sum(.) > 0)) %>%
  left_join(sample_metadata[c(1,10)], by = "Tube_code")
```

```{r commun_wilcox_elem_warm_acc, comment="", message=FALSE, warning=FALSE}
significant_elements_warm <- element_gift_warm %>%
  pivot_longer(-c(Tube_code,time_point), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ time_point)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)%>%
  rownames_to_column(., "Elements")  %>% 
  select(-Elements)

element_gift_sig_warm <- element_gift_warm %>%
  select(Tube_code, all_of(intersect(
    significant_elements_warm$trait,
    colnames(element_gift_warm)
  ))) %>%
  left_join(., sample_metadata[c(1, 10)], by = join_by(Tube_code == Tube_code))

difference_table_warm <- element_gift_sig_warm %>%
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
  mutate(Difference=Acclimation-Wild)%>% 
  mutate(group_color = ifelse(Difference <0, "Wild","Acclimation")) 
```

```{r commun_wilcox_elem_plot_hot, comment="", message=FALSE, warning=FALSE, fig.height=6, fig.width=12, fig.fullwidth=TRUE}
difference_table_warm %>%
  ggplot(aes(x=forcats::fct_reorder(Function,Difference), y=Difference, fill=group_color)) + 
  geom_col() +
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





#NMDS plots (not include in chapter 4)
***NMDS***
  ```{r beta_neutral_nmds_post3_comparison, warning=FALSE, comments="", message=FALSE, results="hide"}
beta_richness_nmds_post3 <- richness_post3 %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace = FALSE) %>%
  scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(subset_meta_post3, by = join_by(sample == Tube_code)) %>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Treatment"),
                     labels=c("Cold-Cold", "Cold-Hot"),
                     values=c("#4477AA","#76b183")) +
  geom_point(size=2) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y = "NMDS2", x="NMDS1 \n Richness beta diversity") +
  theme_classic() +
  theme(legend.position="none")

beta_neutral_nmds_post3 <- neutral_post3 %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace = FALSE) %>%
  scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(subset_meta_post3, by = join_by(sample == Tube_code))%>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Treatment"),
                     labels=c("Cold-Cold", "Cold-Hot"),
                     values=c("#4477AA","#76b183")) +
  geom_point(size=2) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic() +
  theme(legend.position="none")

beta_phylogenetic_nmds_post3 <- phylo_post3 %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace = FALSE) %>%
  scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(subset_meta_post3, by = join_by(sample == Tube_code))%>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Treatment"),
                     labels=c("Cold-Cold", "Cold-Hot"),
                     values=c("#4477AA","#76b183")) +
  geom_point(size=2) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y= element_blank (), x="NMDS1 \n Phylogenetic beta diversity") +
  theme_classic() +
  theme(legend.position="none")
```
