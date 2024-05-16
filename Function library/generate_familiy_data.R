######################################################
#                                                    #
# This code constructs a test family dataset for     #
# testing the fisher scoring algorithms              #
#                                                    #
######################################################



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


# Specify number of individuals in the family
n_individuals <- 10

# Specify id of each individual in the family
id  <- 1:n_individuals

# Randomly assign sex to each individual in the family
sex <- sample(c(1,2), n_individuals, prob = c(1/2,1/2), replace = T)

# Get ids of men and women in the id column
men <- id[sex == 1]
women <- id[sex == 2]

# Pair random men with random women and store the remaining as singles
first_dad <- sample(men, 1)
first_mom <- sample(women,1)

# Remove first parents from potential children
potential_male_children <- men[men != first_dad]
potential_female_children <- women[women != first_mom]

while (potential_male_children > 0 & potential_female_children > 0){
  
  # Assign potential children to parents
  no_male_children <- sample(0:length(potential_male_children), 1)
  no_female_children <- sample(0:length(potential_female_children), 1)
  
  
}



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










