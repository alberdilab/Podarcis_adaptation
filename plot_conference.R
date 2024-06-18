##Plots for the conference poster

#1- wild beta neutral diversity

beta_q1n_nmds_wild %>%
  group_by(Population) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=Population, shape=type)) +
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080',"#d57d2c")) +
  scale_shape_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("circle","square","triangle")) +
  geom_point(size=3) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2, size=0.8) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic()+
  theme(text = element_text(size = 20),
        (legend.text=element_text(size=20))) 
  

#2-fmt vs post fmt neutral beta diversity


#diverging plot

#Create newID to identify duplicated samples
transplants_metadata<-sample_metadata%>%
  mutate(Tube_code=str_remove_all(Tube_code, "_a"))
transplants_metadata$newID <- paste(transplants_metadata$Tube_code, "_", transplants_metadata$individual)

transplant4<-transplants_metadata%>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2" |time_point == "0_Wild")%>%
  column_to_rownames("newID")

transplant4_pairwise <- transplants_metadata %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2"|time_point == "0_Wild")


full_counts<-temp_genome_counts %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Tube_code")%>%
  full_join(transplants_metadata,by = join_by(Tube_code == Tube_code))

transplant4_counts<-full_counts %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2"|time_point == "0_Wild") %>%
  subset(select=-c(315:324)) %>%
  column_to_rownames("newID")%>%
  subset(select=-c(1))%>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)

identical(sort(colnames(transplant4_counts)),sort(as.character(rownames(transplant4))))


beta_div_neutral_transplant4<-hillpair(data=transplant4_counts, q=1)

transplant4_beta_div_neutral<-as.matrix(beta_div_neutral_transplant4$S) #convert the pairwise comparisons into a matrix

transplant4_beta_div_neutral_df <- as.data.frame(as.table(as.matrix(transplant4_beta_div_neutral))) #convert the matrix to a data frame with pairwise values

transplant4_beta_div_neutral_df <- transplant4_beta_div_neutral_df %>% #remove the duplicated pairwise comparisons
  filter(Var1 != Var2) %>%
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(..1, ..2)), collapse = "_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(Sample_1 = Var1, Sample_2 = Var2, Distance = Freq) %>%
  arrange(Sample_1, Sample_2)


transplant4_beta_div_neutral_df <- transplant4_beta_div_neutral_df %>% #keep only the pairwise comparisons between the same individual
  mutate(Sample_1_part = sub("^[^_]*_", "", Sample_1),
         Sample_2_part = sub("^[^_]*_", "", Sample_2)) %>%
  filter(Sample_1_part == Sample_2_part) %>%
  select(Sample_1, Sample_2, Distance)


transplant4_beta_div_neutral_met<-transplant4_beta_div_neutral_df %>% #merge with the metadata associated to the pairwise comparisons for Sample_1
  inner_join(transplant4_pairwise, by = join_by(Sample_1 == "newID"))

# Extract the part before the first `_` for Sample_1 and Sample_2
transplant4_beta_div_neutral_met$Sample_1_part <- sub("_.*", "", transplant4_beta_div_neutral_met$Sample_1)

# Create a named vector for matching Tube_code with Time_point
time_point_map <- setNames(transplant4_beta_div_neutral_met$time_point, transplant4_beta_div_neutral_met$Tube_code)

# Function to replace the part before the first _ with the corresponding Time_point
replace_with_time_point <- function(sample_part, tube_code, time_point_map) {
  if (tube_code %in% names(time_point_map)) {
    return(time_point_map[tube_code])
  } else {
    return(sample_part)
  }
}

# Apply the function to Sample_1_part
transplant4_beta_div_neutral_met$Sample_1_part <- mapply(replace_with_time_point,
                                                         transplant4_beta_div_neutral_met$Sample_1_part,
                                                         transplant4_beta_div_neutral_met$Tube_code,
                                                         list(time_point_map))

transplant4_beta_div_neutral_met<-transplant4_beta_div_neutral_met %>%
  subset(select=-c(6:16))

transplant4_beta_div_neutral_met$Sample_2_part <- sub("_.*", "", transplant4_beta_div_neutral_met$Sample_2)


transplant4_beta_div_neutral_met<-transplant4_beta_div_neutral_met %>% #merge with the metadata associated to the pairwise comparisons for Sample_2
  inner_join(transplant4_pairwise, by = join_by(Sample_2 == "newID"))

# Create a named vector again for matching Tube_code with Time_point
time_point_map <- setNames(transplant4_beta_div_neutral_met$time_point, transplant4_beta_div_neutral_met$Tube_code)


replace_with_time_point <- function(sample_part, tube_code, time_point_map) {
  if (tube_code %in% names(time_point_map)) {
    return(time_point_map[tube_code])
  } else {
    return(sample_part)
  }
}

# Apply the function to Sample_2_part
transplant4_beta_div_neutral_met$Sample_2_part <- mapply(replace_with_time_point,
                                                         transplant4_beta_div_neutral_met$Sample_2_part,
                                                         transplant4_beta_div_neutral_met$Tube_code,
                                                         list(time_point_map))


# Create a new variable to plot the distances between the time_points
transplant4_beta_div_neutral_met$comparison <- paste(transplant4_beta_div_neutral_met$Sample_1_part, "-",transplant4_beta_div_neutral_met$Sample_2_part)


##plot the differences between wild and post-transplant and transplant and post-transplant

transplant4_beta_div_neutral_met  %>% 
  ggplot(aes(x=type, y=Distance, color=type, fill=type))+
    geom_boxplot() +
    scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
    scale_fill_manual(name="Type",
                      breaks=c("Control", "Hot_control", "Treatment"),
                      labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                      values=c("#4477AA50","#d57d2c50","#76b18350")) +
    scale_x_discrete(labels=c("Control" = "Cold-Cold", "Hot_control" = "Hot-Hot", "Treatment" = "Cold-Hot")) +
    theme_minimal() +
    facet_wrap(~comparison, labeller = labeller(comparison = c("0_Wild - 6_Post-FMT2" = "Wild-Treatment", "0_Wild - 4_Transplant2" = "Wild-Transplant", "6_Post-FMT2 - 4_Transplant2"="Treatment-Transplant"))) +
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1))



#individual comparison

transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  mutate(total_dis=sum(Distance)) %>%
  filter(sample_n == 2) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(rel_dis=Distance/total_dis) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=individual,x=rel_dis,fill=comparison))+
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(name="Timeline comparison",
                    breaks=c("0_Wild - 6_Post-FMT2","6_Post-FMT2 - 4_Transplant2"),
                    labels=c("Wild vs Treatment","Treatment vs Transplant"),
                    values=c("#e38888","#7B9381" )) +
  facet_wrap(.~ type, scales="free", ncol=1)+
  labs(y = "Individual", x="Real distance")+
  theme(
    axis.text.y = element_text("none"))
  


#try to make another type of barplot
transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  mutate(total_dis=sum(Distance)) %>%
  filter(sample_n == 2) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=type,x=Distance, type))+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(name="Timeline comparison",
                    breaks=c("0_Wild - 6_Post-FMT2","6_Post-FMT2 - 4_Transplant2"),
                    labels=c("Wild vs Treatment","Treatment vs Transplant"),
                    values=c("#e38888","#7B9381" )) +
  facet_wrap(.~ comparison, scales="free", ncol=1)+
  labs(y = "Individual", x="Real distance")



h<-transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  filter(sample_n == 2) %>%
  mutate(total_dis=sum(Distance)) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=comparison,x=Distance,fill=type))+
  geom_col(position=position_dodge())  +
  facet_wrap(.~ type, ncol=1)+
  labs(y = "Time-point comparison", x="Pairwise dissimilarity distance")



#stripplot
transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  mutate(total_dis=sum(Distance)) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=Distance,x=type,fill=type,))+
  geom_jitter(
  aes(shape = type, color = type), size = 1.2,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  stat_summary(
    aes(color = type), fun.data="mean_sdl", fun.args = list(mult=1), 
    size = 0.4, position = position_dodge(0.8)
  )+
  stat_summary(aes(color = type), size = 0.4,
               fun.data="mean_sdl",  fun.args = list(mult=1))+
  facet_wrap(.~ comparison, ncol=1)+
  labs(y = "Group type", x="Pairwise dissimilarity distance")



#Combination of plots

p<-transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  mutate(total_dis=sum(Distance)) %>%
  filter(sample_n == 2) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(rel_dis=Distance/total_dis) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=type,x=rel_dis,fill=comparison))+
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(name="Timeline comparison",
                    breaks=c("0_Wild - 6_Post-FMT2","6_Post-FMT2 - 4_Transplant2"),
                    labels=c("Wild vs Treatment","Treatment vs Transplant"),
                    values=c("#e38888","#7B9381" )) +
  facet_wrap(.~ type, scales="free", ncol=1)+
  labs(y = "Individual", x="Real distance")+
  theme(
    axis.text.y = element_text("none"), legend.position = "bottom")
  
q<-transplant4_beta_div_neutral_met  %>% 
  filter(comparison=="0_Wild - 4_Transplant2")  %>% 
  ggplot(aes(x=type, y=Distance, color=type, fill=type))+
  geom_boxplot() +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("#4477AA","#d57d2c","#76b183")) +
  scale_fill_manual(name="Type",
                    breaks=c("Control", "Hot_control", "Treatment"),
                    labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                    values=c("#4477AA50","#d57d2c50","#76b18350")) +
  scale_x_discrete(labels=c("Control" = "Cold-Cold", "Hot_control" = "Hot-Hot", "Treatment" = "Cold-Hot")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  labs(title="Wild vs Transplant",y = "Pairwise distance", x="Type")
  

ggarrange(q, p, ncol=2, nrow=1, legend="right")



#plot combination 2

h<-transplant4_beta_div_neutral_met %>%
  filter(comparison %in% c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2")) %>%
  group_by(individual) %>%
  mutate(sample_n=n()) %>%
  filter(sample_n == 2) %>%
  mutate(total_dis=sum(Distance)) %>%
  select(comparison,Distance,total_dis,individual,type) %>%
  arrange(individual) %>%
  mutate(comparison=factor(comparison,levels=c("6_Post-FMT2 - 4_Transplant2","0_Wild - 6_Post-FMT2"))) %>%
  ggplot(aes(y=comparison,x=Distance,fill=type))+
  geom_col(position=position_dodge())  +
  facet_wrap(.~ type, ncol=1)+
  labs(y = "Time-point comparison", x="Pairwise dissimilarity distance")

ggarrange(q, h, ncol=2, nrow=1, legend="right")


#transplant vs post transplant comparison
beta_neutral_nmds_transplant3 %>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("#76b183","#d57d2c", "#4477AA"))+
  scale_shape_manual(name="Time-point",
                     breaks=c("4_Transplant2", "6_Post-FMT2"),
                     labels=c("Transplant", "Post-Transplant"),
                     values=c(16, 17))+
  geom_point(size=2) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic()


beta_neutral_nmds_transplant3 %>%
  group_by(individual) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("#76b183","#d57d2c", "#4477AA")) +
  scale_shape_manual(name="Time-point",
                     breaks=c("4_Transplant2", "6_Post-FMT2"),
                     labels=c("Transplant", "Post-Transplant"),
                     values=c(16, 17)) +
  geom_point(size=2) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic()

#plotting as ellipse
beta_neutral_nmds_transplant3 %>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("#76b183","#d57d2c", "#4477AA")) +
  scale_shape_manual(name="Time-point",
                     breaks=c("4_Transplant2", "6_Post-FMT2"),
                     labels=c("Transplant", "Post-Transplant"),
                     values=c(16, 17)) +
  #stat_ellipse()+
  geom_point(size=2) +
  #geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic()


##adonis comparison of each type between the time_points

###control

transplant3_control<-transplants_metadata%>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2")%>%
  filter(type=="Control")%>%
  column_to_rownames("newID")

transplant3_counts_control<-full_counts %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2") %>%
  filter(type=="Control") %>%
  subset(select=-c(315:324)) %>%
  column_to_rownames("newID")%>%
  subset(select=-c(1))%>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)

beta_div_neutral_transplant3_control<-hillpair(data=transplant3_counts_control, q=1)

adonis2(formula=beta_div_neutral_transplant3_control$S ~ time_point, data=transplant3_control[labels(beta_div_neutral_transplant3_control$S),], permutations=999) %>%
  as.matrix() %>%
  kable()

pairwise<-pairwise.adonis(beta_div_neutral_transplant3_control$S,transplant3_control$time_point, perm=999)
pairwise

###hot
transplant3_hot<-transplants_metadata%>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2")%>%
  filter(type=="Hot_control")%>%
  column_to_rownames("newID")

transplant3_counts_hot<-full_counts %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2") %>%
  filter(type=="Hot_control") %>%
  subset(select=-c(315:324)) %>%
  column_to_rownames("newID")%>%
  subset(select=-c(1))%>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)

beta_div_neutral_transplant3_hot<-hillpair(data=transplant3_counts_hot, q=1)

adonis2(formula=beta_div_neutral_transplant3_hot$S ~ time_point, data=transplant3_hot[labels(beta_div_neutral_transplant3_hot$S),], permutations=999) %>%
  as.matrix() %>%
  kable()

pairwise<-pairwise.adonis(beta_div_neutral_transplant3_hot$S,transplant3_hot$time_point, perm=999)
pairwise

###treatment
transplant3_treatment<-transplants_metadata%>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2")%>%
  filter(type=="Treatment")%>%
  column_to_rownames("newID")

transplant3_counts_treatment<-full_counts %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2") %>%
  filter(type=="Treatment") %>%
  subset(select=-c(315:324)) %>%
  column_to_rownames("newID")%>%
  subset(select=-c(1))%>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)

beta_div_neutral_transplant3_treatment<-hillpair(data=transplant3_counts_treatment, q=1)

adonis2(formula=beta_div_neutral_transplant3_treatment$S ~ time_point, data=transplant3_treatment[labels(beta_div_neutral_transplant3_treatment$S),], permutations=999) %>%
  as.matrix() %>%
  kable()

pairwise<-pairwise.adonis(beta_div_neutral_transplant3_treatment$S,transplant3_treatment$time_point, perm=999)
pairwise

#3-post fmt neutral beta diversity

beta_neutral_nmds_post2 %>%
  group_by(type) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
  scale_color_manual(name="Type",
                     breaks=c("Control", "Hot_control", "Treatment"),
                     labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                     values=c("#4477AA","#d57d2c","#76b183")) +
  geom_point(size=3) +
  geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2, size=0.8) +
  labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
  theme_classic()+
  theme(text = element_text(size = 20),
        (legend.text=element_text(size=20))) 



##Plot PhD day

GIFTs_functions_community %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  select(c(1:21, 24)) %>%
  pivot_longer(-c(sample,Population),names_to = "trait", values_to = "value") %>%
  mutate(trait = case_when(
    trait %in% GIFT_db$Code_function ~ GIFT_db$Function[match(trait, GIFT_db$Code_function)],
    TRUE ~ trait
  )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=value, y=Population, group=Population, fill=Population, color=Population)) +
  geom_boxplot() +
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080',"#d57d2c")) +
  scale_fill_manual(name="Population",
                    breaks=c("Cold_wet","Hot_dry"),
                    labels=c("Cold","Hot"),
                    values=c('#00808050',"#d57d2c50")) +
  scale_y_discrete(breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot")) +
  facet_grid(trait ~ ., space="free", scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(angle = 0)) + 
  labs(y="Traits",x="Metabolic capacity index")
