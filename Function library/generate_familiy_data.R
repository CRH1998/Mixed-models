######################################################
#                                                    #
# This code constructs a test family dataset for     #
# testing the fisher scoring algorithms              #
#                                                    #
######################################################

#Construct function that pairs parents with children
pair_parents_with_children <- function(parents, unassigned_individuals){
  
  sub_family_df <- data.frame()
  
  for (i in 1:nrow(parents)){
    # Assign potential children to parents
    no_children <- rbinom(1,nrow(unassigned_individuals),0.25)
    children <- unassigned_individuals[sample(nrow(unassigned_individuals), no_children),]
    
    # Construct family dataframe
    ifelse(no_children == 0, sub_family_df <- rbind(sub_family_df, data.frame(NA, NA, NA, NA)), 
           sub_family_df <- rbind(sub_family_df, data.frame(children, parents[i,1], parents[i,3])))
    
    #Remove children who are already assigned
    unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% children$id),]
  }
  names(sub_family_df) <- c("id", "sex","dadid", "momid")
  return(sub_family_df)
}


# Construct function that pairs children with spouses
pair_children_with_spouses <- function(family_df, unassigned_individuals){
  
  # Get sex of children
  male_children <- family_df[family_df$sex == 1,]
  female_children <- family_df[family_df$sex == 2,]
  
  # Get ids of men and women in the id column
  unassigned_men <- unassigned_individuals[unassigned_individuals$sex == 1,]
  unassigned_women <- unassigned_individuals[unassigned_individuals$sex == 2,]
  
  # Sample spouses for children among the remaining
  no_male_spouses <- sample(0:min(nrow(female_children), nrow(unassigned_men)), 1)
  no_female_spouses <- sample(0:min(nrow(male_children), nrow(unassigned_women)), 1)
  
  married_male_unassigned <- unassigned_men[sample(nrow(unassigned_men), no_male_spouses),]
  married_female_unassigned <- unassigned_women[sample(nrow(unassigned_women), no_female_spouses),]
  
  # Sample children:
  married_male_children <- male_children[sample(nrow(male_children), no_female_spouses),][,c("id", "sex")]
  married_female_children <- female_children[sample(nrow(female_children), no_male_spouses),][,c("id", "sex")]
  
  # Construct parent dataframe
  parents <- data.frame(rbind(cbind(married_male_children, married_female_unassigned), 
                              cbind(married_male_unassigned, married_female_children)))
  
  if (no_male_spouses == 0 & no_female_spouses == 0){
    assigned_individuals <- data.frame(matrix(ncol = 4, nrow = 0))
    names(assigned_individuals) <- c("id", "sex","dadid", "momid")
    return(list('parents' = parents, 'assigned_individuals' = assigned_individuals))
  }
  
  assigned_individuals <- data.frame(rbind(married_male_unassigned, married_female_unassigned), 0, 0)
  names(assigned_individuals) <- c("id", "sex","dadid", "momid")
  
  return(list('parents' = parents, 'assigned_individuals' = assigned_individuals))
}

set.seed(1)

# Specify number of individuals in the family
n_individuals <- 10

# Specify id of each individual in the family
id  <- 1:n_individuals

# Randomly assign sex to each individual in the family
sex <- sample(c(1,2), n_individuals, prob = c(1/2,1/2), replace = T)

# Unassigned individuals
unassigned_individuals <- data.frame(id,sex)

# Get ids of men and women in the id column
men <- unassigned_individuals[unassigned_individuals$sex == 1,]
women <- unassigned_individuals[unassigned_individuals$sex == 2,]

# Pair random man with random woman to create first parents
first_dad <- men[sample(nrow(men), 1),]
first_mom <- women[sample(nrow(women), 1),]

parents <- data.frame(first_dad, first_mom)

# Remove first parents from unassigned individuals
unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% c(first_dad$id, first_mom$id)),]

final_family_df <- data.frame(rbind(first_dad, first_mom), 0, 0)
names(final_family_df) <- c("id", "sex","dadid", "momid")

while(nrow(unassigned_individuals) > 0){
  
  set.seed(1)
  
  if(nrow(parents) > 0){
    family_df <- pair_parents_with_children(parents, unassigned_individuals)
  }
  # Update unassigned individuals
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% family_df$id),]
  
  if (((sum(unassigned_individuals$sex == 1) == 0) | (sum(unassigned_individuals$sex == 2) == 0)) & (nrow(unassigned_individuals) > 0)){
    last_unassinged <- cbind(unassigned_individuals, family_df[nrow(family_df),3], family_df[nrow(family_df),4])
    names(last_unassinged) <- c("id", "sex","dadid", "momid")
    final_family_df <- rbind(final_family_df, family_df, last_unassinged)
    break
  } else if (nrow(parents) > 0) {
    final_family_df <- rbind(final_family_df, family_df)
  }
  
  if (nrow(unassigned_individuals) > 0  & nrow(family_df) > 0){
    parents <- pair_children_with_spouses(family_df, unassigned_individuals)$parents
    assigned_individuals <- pair_children_with_spouses(family_df, unassigned_individuals)$assigned_individuals
    final_family_df <- rbind(final_family_df, assigned_individuals)
  }
  
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% assigned_individuals$id),]
}
final_family_df









family_structure_df <- data.frame(id = c(1,2,3,4,5,6,7,8,9,10), 
                                  sex = c(1,2,2,1,1,2,2,2,2,2), 
                                  dadid = c(0,0,1,1,1,0,0,4,4,5), 
                                  momid = c(0,0,2,2,2,0,0,6,6,7), 
                                  famid = 1)
family_structure_df


relation_matrix <- matrix(c(1,2,4,1, 4,6,4,1, 5,7,4,1), ncol = 4, byrow = T)
relation_matrix

foo <- pedigree(id = family_structure_df$id, 
                dadid = family_structure_df$dadid, 
                momid = family_structure_df$momid, 
                sex = family_structure_df$sex,
                famid = family_structure_df$famid,
                relation = relation_matrix)
ped <- foo['1']
plot(ped)








while (potential_male_children > 0 & potential_female_children > 0){
  
  # Assign potential children to parents
  no_male_children <- sample(0:length(potential_male_children), 1)
  no_female_children <- sample(0:length(potential_female_children), 1)
  
  
}

















relation_matrix <- matrix(c(1,2,4,1, 4,6,4,1, 5,7,4,1), ncol = 4, byrow = T)


dadid <- sample(c(men,0), size = n_individuals, replace = T)
momid <- sample(c(women,0), size = n_individuals, replace = T)

dadid[dadid == id] <- 0
momid[momid == id] <- 0

famid <- 1

pedigree(id = id, 
         dadid = dadid, 
         momid = momid, 
         sex = sex,
         famid = famid)









male_children <- sample(potential_male_children, no_male_children)
female_children <- sample(potential_female_children, no_female_children)


#Constructing first part of family dataframe
rbind(cbind(male_children, first_dad, first_mom), cbind(female_children, first_dad, first_mom))


remaining_men <- potential_male_children[potential_male_children != male_children]
remaining_women <- potential_female_children[potential_female_children != female_children]

# Sampling spouses for children among the remaining men and women
sampling_male_spouses <- sample(remaining_men, sample(0:length(female_children),1))
sampling_female_spouses <- sample(remaining_women, sample(0:length(male_children),1))

# Coupling children with spouses
second_dads <- sample(male_children, length(sampling_female_spouses))
second_moms <- sample(female_children, length(sampling_male_spouses))




