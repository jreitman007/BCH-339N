get_node_value <- function(node1, node2) {
  differences <- 0
  
  for (i in 1:(length(node1))) {
    s1 <- sort(unlist(strsplit(node1[i], ",")))
    s2 <- sort(unlist(strsplit(node2[i], ",")))
    
    if (length(intersect(s1,s2)) == 0) {
      differences <- differences + 1
    }
  }
  probability <-  exp(-8*(differences/length(node1)))
  
  return(probability)
}


get_consensus_sequence <- function(sequence1, sequence2) {
  new_sequence <- c()
  
  for (i in 1:(length(sequence1))) {
    s1 <- sort(unlist(strsplit(sequence1[i], ",")))
    s2 <- sort(unlist(strsplit(sequence2[i], ",")))
    
    if (length(intersect(s1,s2)) > 0) {
      new_sequence[i] <- paste(intersect(s1, s2), collapse = ",")
    } else {
      new_sequence[i] <- paste(union(s1,s2), collapse = ",")
    }
  }
  
  return(new_sequence)
}

tree_likelihood <- function(sequences) {
  probability <- c()
  new_probability <- c()
  new_sequences <- c()
  
  number <- length(sequences)  
  
  for (i in  2*1:(length(sequences)/2)) {
    probability[i/2] <- get_node_value(sequences[[i]], sequences[[i-1]])
    
    new_sequences[[i/2]] <- get_consensus_sequence(sequences[[i]], sequences[[i - 1 ]])
  }
  
  if (number %% 2 == 1) {
    probability[i / 2] <- probability[i / 2] * get_node_value(new_sequences[[i / 2]], sequences[number])
  }
  
  sequences <- new_sequences
  
  while(length(probability) != 1) {
    number <- length(sequences)  
    new_sequences <- c()
    new_probability <- c()
    for (i in 2*1:(length(probability)/2)) {
      
      new_probability[i/2] <- probability[i] * probability[i - 1] * get_node_value(sequences[[i]], sequences[[i - 1]])
      
      new_sequences[[i/2]] <- get_consensus_sequence(sequences[[i]], sequences[[i - 1]])
    }
    
    if (number %% 2 == 1) {
      new_probability[i / 2] <- new_probability[i / 2] * probability[number] * get_node_value(new_sequences[[i / 2]], sequences[[number]])
      new_sequences[[i/2]] <- get_consensus_sequence(new_sequences[[i / 2]], sequences[[number]])  
    }
    probability <- new_probability
    
    sequences <- new_sequences
    
  } 
  
  return(probability)
}

generate_trees <- function(sequences, names) {
  
  probs <- 0
  
  candidate_sequences <- sequences
  candidate_names <- names
  best_names <- names
  best_sequence <- sequences
  
  for (i in 1:length(sequences)) {
    for (j in 1:length(sequences)) {
      candidate_names <- best_names
      candidate_sequences <- best_sequence
      
      if( i != j) {
        
        print("Best")
        print(best_names)
        sequence1 <- candidate_sequences[[i]]
        sequence2 <- candidate_sequences[[j]]
        name1 <- candidate_names[i]
        name2 <- candidate_names[j]
        
        
        candidate_sequences[[i]] <- sequence2
        candidate_sequences[[j]] <- sequence1
        candidate_names[i] <- name2
        candidate_names[j] <- name1
        
        new_prob <- tree_likelihood(candidate_sequences)
        
        print("Candidate")
        print(candidate_names)
        print("Likelihood")
        print(toString(new_prob))
        
        if (new_prob > probs) {
          probs <- new_prob
          print(toString(probs))
          
          best_sequence <- candidate_sequences
          best_names <- candidate_names
          
          i <- 1
          j <- 1
        }
        
      }
      
    }  
  }
  
  print("Best Likelihood")
  print(toString(probs))
  print("Tree")
  print(best_names)
  return(best_names)
}



#Test data

sequences <- list(DQ182595, JX869059, KT368829, GU553363, DQ648857, DQ412043, MT163719, MT163718, MT135043, MT135041)
names <- c("DQ182595", "JX869059", "KT368829", "GU553363", "DQ648857", "DQ412043",  "MT163719", "MT163718", "MT135043", "MT135041")

for (i in 1:10) {
  sequences[[i]] <- replace(sequences[[i]], sequences[[i]] == "-", "A,C,G,T")
}

start_time <- Sys.time()
generate_trees(sequences, names)
end_time <- Sys.time()