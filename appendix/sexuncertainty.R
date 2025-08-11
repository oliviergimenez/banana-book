

####################################################  sex uncertainty

# STATES
# alive as male (M)
# alive as female (F)
# dead (D)

# OBS
# not seen (1)
# judged from copulation to be M (2)
# judged from begging food to be M (3)
# judged from coutship feeding to be M (4)
# judged from body size to be M (5)
# judged from copulation to be F (6)
# judged from begging food to be F (7)
# judged from coutship feeding to be F (8)
# judged from body size to be F (9)
# not judged (10)

# # Read the data (as a matrix, assuming space-separated values in a file)
# dat <- as.matrix(read.table("/Users/oliviergimenez/Dropbox/OG/GITHUB/book/matos/case-studies/covariates/sex-uncertainty/sexo_audouines9/Sexaudouin9.text"))  # or use here::here("path", "to", "file.txt")
# 
# # Separate sequence and count
# sequence <- dat[, 1:(ncol(dat)-1)]
# counts <- dat[, ncol(dat)]
# 
# ind <- c(1, 176, 991, 177, 178, 179, 468, 667, 668, 669, 180, 670, 469, 671, 181, 182, 830, 183, 917, 184, 918, 992, 672, 993) 
# length(ind)
# 
# 
# 
# geneticsex <- c(rep('F', 11), rep('M', 13))
# 
# # Expand rows according to counts
# expanded <- sequence[rep(1:nrow(sequence), counts), ]
# 
# # Convert to data.frame (optional)
# expanded_df <- as.data.frame(expanded)
# 
# # Write to file (optional)
# write.table(expanded_df, file = "audouin-gull.txt", row.names = FALSE, col.names = FALSE)
# 

# read in data
y <- as.matrix(read.table(here::here("dat", "audouin-gull.txt")))

library(nimble)

hmm.sex <- nimbleCode({
  
  # priors
  phiF ~ dunif(0, 1) # prior survival male
  phiM ~ dunif(0, 1) # prior survival female
  p ~ dunif(0, 1) # prior encounter prob
  pi ~ dunif(0, 1) # prob init state male
  e ~ dunif(0,1)
  for (j in 1:4){
    x[j] ~ dunif(0,1)
  }
  m[1] ~ dunif(0,1)
  m[2] <- (1 - m[1]) * alpha
  alpha ~ dunif(0,1)
  m[3] <- 1 - m[1] - m[2]
  m[4] ~ dunif(0,1)
  
  # HMM ingredients
  delta[1] <- pi
  delta[2] <- 1 - pi
  delta[3] <- 0
  
  gamma[1,1] <- phiM            # Pr(C t -> C t+1)
  gamma[1,2] <- 0               # Pr(C t -> A t+1)
  gamma[1,3] <- 1 - phiM        # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0               # Pr(A t -> C t+1)
  gamma[2,2] <- phiF            # Pr(A t -> A t+1)
  gamma[2,3] <- 1 - phiF        # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0               # Pr(dead t -> alive t+1)
  gamma[3,2] <- 0               # Pr(dead t -> alive t+1)
  gamma[3,3] <- 1               # Pr(dead t -> dead t+1)
  
  omega[1,1] <- 1 - p                                  
  omega[1,2] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[1,3] <- p * e * (1 - m[4]) * m[2] * x[2]      
  omega[1,4] <- p * e * (1 - m[4]) * m[3] * x[3]      
  omega[1,5] <- p * e * m[4] * x[4]      
  omega[1,6] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[1,7] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[1,8] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[1,9] <- p * e * m[4] * (1 - x[4])      
  omega[1,10] <- p * (1 - e)      
  omega[2,1] <- 1 - p                                  
  omega[2,2] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[2,3] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[2,4] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[2,5] <- p * e * m[4] * (1 - x[4])      
  omega[2,6] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[2,7] <- p * e * (1 - m[4]) * m[2] * x[2]     
  omega[2,8] <- p * e * (1 - m[4]) * m[3] * x[3]   
  omega[2,9] <- p * e * m[4] * x[4]   
  omega[2,10] <- p * (1 - e)      
  omega[3,1] <- 1                                  
  omega[3,2] <- 0     
  omega[3,3] <- 0     
  omega[3,4] <- 0     
  omega[3,5] <- 0     
  omega[3,6] <- 0     
  omega[3,7] <- 0     
  omega[3,8] <- 0     
  omega[3,9] <- 0     
  omega[3,10] <- 0     
  
  omega.init[1,1] <- 0                                 
  omega.init[1,2] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[1,3] <- e * (1 - m[4]) * m[2] * x[2]      
  omega.init[1,4] <- e * (1 - m[4]) * m[3] * x[3]      
  omega.init[1,5] <- e * m[4] * x[4]      
  omega.init[1,6] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[1,7] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[1,8] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[1,9] <- e * m[4] * (1 - x[4])      
  omega.init[1,10] <- (1 - e)      
  omega.init[2,1] <- 0                                  
  omega.init[2,2] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[2,3] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[2,4] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[2,5] <- e * m[4] * (1 - x[4])      
  omega.init[2,6] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[2,7] <- e * (1 - m[4]) * m[2] * x[2]     
  omega.init[2,8] <- e * (1 - m[4]) * m[3] * x[3]   
  omega.init[2,9] <- e * m[4] * x[4]   
  omega.init[2,10] <- (1 - e)      
  omega.init[3,1] <- 1                                  
  omega.init[3,2] <- 0     
  omega.init[3,3] <- 0     
  omega.init[3,4] <- 0     
  omega.init[3,5] <- 0     
  omega.init[3,6] <- 0     
  omega.init[3,7] <- 0     
  omega.init[3,8] <- 0     
  omega.init[3,9] <- 0     
  omega.init[3,10] <- 0     
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:10])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:10])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first) 
my.data <- list(y = y + 1)


# # Function to create initial latent states z
# #create_initial_z <- function(y, first) {
#   N <- nrow(y)
#   K <- ncol(y)
#   zinit <- matrix(NA, nrow = N, ncol = K)
#   
#   for (i in 1:N) {
#     # Détermination du sexe à la première détection
#     first_obs <- y[i, first[i]]  # valeurs déjà de 0-9, pas besoin de +1
#     
#     # Classification basée sur les codes d'observation
#     # 0 = non observé, 1-5 = mâle, 6-8 = femelle, 9 = incertain
#     if (first_obs %in% 1:5) {
#       initial_sex <- 1  # mâle
#     } else if (first_obs %in% 6:8) {
#       initial_sex <- 2  # femelle
#     } else {
#       initial_sex <- sample(1:2, 1)  # sexe aléatoire si incertain (9)
#     }
#     
#     zinit[i, first[i]] <- initial_sex
#     
#     # Pour les occasions suivantes
#     current_state <- initial_sex
#     if (first[i] < K) {
#       for (j in (first[i] + 1):K) {
#         obs <- y[i, j]
#         
#         if (obs == 0) {  # non observé
#           # Probabilité de survie élevée si pas d'évidence de mort
#           if (runif(1) < 0.85) {
#             zinit[i, j] <- current_state  # survit
#           } else {
#             zinit[i, j] <- 3  # meurt
#             current_state <- 3
#           }
#         } else {  # observé vivant (1-9)
#           # Maintient l'état vivant
#           if (current_state == 3) {
#             # Si était mort, ne peut pas revenir vivant - problème dans les données
#             current_state <- sample(1:2, 1)
#           }
#           zinit[i, j] <- current_state
#         }
#       }
#     }
#   }
#   
#   return(zinit)
# }

# # Function to create initial latent states z
create_initial_z <- function(y, first) {
  N <- nrow(y)
  K <- ncol(y)
  zinit <- matrix(NA, nrow = N, ncol = K)
  
  for (i in 1:N) {
    # Détermination du sexe à la première détection
    first_obs <- y[i, first[i]]  # valeurs déjà de 0-9, pas besoin de +1
    
    # Classification basée sur les codes d'observation
    # 0 = non observé, 1-5 = mâle, 6-8 = femelle, 9 = incertain
    if (first_obs %in% 1:5) {
      initial_sex <- 1  # mâle
    } else if (first_obs %in% 6:8) {
      initial_sex <- 2  # femelle
    } else {
      initial_sex <- sample(1:2, 1)  # sexe aléatoire si incertain (9)
    }
    
    zinit[i, first[i]] <- initial_sex
    
    # Pour les occasions suivantes
    current_state <- initial_sex
    if (first[i] < K) {
      for (j in (first[i] + 1):K) {
        obs <- y[i, j]
        
        if (obs == 0) {  # non observé
          zinit[i, j] <- current_state  # survit
        } else {  # observé vivant (1-9)
          # Maintient l'état vivant
          if (current_state == 3) {
            # Si était mort, ne peut pas revenir vivant - problème dans les données
            current_state <- sample(1:2, 1)
          }
          zinit[i, j] <- current_state
        }
      }
    }
  }
  
  return(zinit)
}

# Function to generate valid m parameters
generate_valid_m <- function() {
  m <- rep(NA, 4)
  
  # m[4] : proportion de la 4ème catégorie
  m[4] <- runif(1, 0.1, 0.4)
  
  # m[1] et m[2] avec contrainte que m[1] + m[2] + m[3] = 1 et m[3] > 0
  repeat {
    m[1] <- runif(1, 0.1, 0.4)
    alpha <- runif(1, 0.1, 0.4)
    m[2] <- (1 - m[1]) * alpha
    if (m[1] + m[2] <= 0.8) break  # assure que m[3] >= 0.2
  }
  m[3] <- 1 - m[1] - m[2]
  
  return(m)
}

# Function to create valid initial parameter values - VERSION CORRIGEE
create_initial_params <- function() {
  # Paramètres avec des valeurs plus réalistes
  phiM <- runif(1, 0.75, 0.95)  # survie mâle
  phiF <- runif(1, 0.75, 0.95)  # survie femelle
  p <- runif(1, 0.4, 0.8)       # probabilité de détection
  pi <- runif(1, 0.45, 0.55)    # ratio des sexes
  e <- runif(1, 0.3, 0.7)       # probabilité d'encounter given detected
  x <- 1 - runif(4, 0.01, 0.3)       # probabilités de sexage correct
  m <- generate_valid_m()
  
  return(list(
    phiM = phiM,
    phiF = phiF, 
    p = p,
    pi = pi,
    e = e,
    x = x,
    m = m
  ))
}

# Create the complete initial values function
initial.values <- function() {
  params <- create_initial_params()
  zinit <- create_initial_z(y, first)
  return(c(params, list(z = zinit)))
}

# # Fonction de diagnostic pour vérifier les valeurs initiales
# check_initial_values <- function() {
#   init_vals <- initial.values()
#   
#   cat("Vérification des paramètres initiaux:\n")
#   cat("phiM:", init_vals$phiM, "\n")
#   cat("phiF:", init_vals$phiF, "\n")
#   cat("p:", init_vals$p, "\n")
#   cat("pi:", init_vals$pi, "\n")
#   cat("e:", init_vals$e, "\n")
#   cat("x:", init_vals$x, "\n")
#   cat("m:", init_vals$m, "\n")
#   cat("Somme m[1:3]:", sum(init_vals$m[1:3]), "\n")
#   
#   # Vérifier quelques probabilités omega
#   e_val <- init_vals$e
#   m_val <- init_vals$m
#   x_val <- init_vals$x
#   p_val <- init_vals$p
#   
#   # Probabilité pour état 1, observation 2
#   prob_example <- p_val * e_val * (1 - m_val[4]) * m_val[1] * x_val[1]
#   cat("Exemple probabilité omega[1,2]:", prob_example, "\n")
#   
#   return(init_vals)
# }

parameters.to.save <- c("phiM", 
                        "phiF", 
                        "p", 
                        "pi",
                        "e", 
                        "m", 
                        "x")

n.iter <- 2500
n.burnin <- 500
n.chains <- 2

# # Diagnostic avant création du modèle
# cat("=== DIAGNOSTIC DES VALEURS INITIALES ===\n")
# test_init <- check_initial_values()
# 
# # Création du modèle avec gestion d'erreur
# cat("\n=== CRÉATION DU MODÈLE ===\n")
# set.seed(123)
# 

# survival <- nimbleModel(code = hmm.sex,
#                         data = my.data,
#                         constants = my.constants,
#                         inits = initial.values())
#   
# survival$calculate()
# 
# for (i in 1:N) {
#   for (j in first[i]:K) {
#     z_ij <- survival$z[i,j]
#     y_ij <- survival$y[i,j]
#     if (z_ij %in% 1:3 && y_ij %in% 1:10) {
#       prob <- survival$omega[z_ij, y_ij]
#       if (prob == 0) {
#         cat("⚠️ omega[", z_ij, ",", y_ij, "] == 0 at y[", i, ",", j, "]\n")
#       }
#     }
#   }
# }
# 
# survival$getNodeNames()
# 
# survival$e
# survival$p
# survival$m
# survival$x
# 
# apply(survival$omega,1,sum)
# apply(survival$omega.init,1,sum)
# 
# survival$gamma
# apply(survival$gamma, 1, sum)
# 
# survival$delta

out <- nimbleMCMC(code = hmm.sex, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

save(out, file = "sex.RData")

library(MCMCvis)

MCMCsummary(out, round = 2)
MCMCtrace(out, pdf = FALSE, params = c("p"))
MCMCtrace(out, pdf = FALSE, params = c("e"))
MCMCtrace(out, pdf = FALSE, params = c("m"))
MCMCtrace(out, pdf = FALSE, params = c("x"))
MCMCtrace(out, pdf = FALSE, params = c("phiM", "phiF"))
MCMCtrace(out, pdf = FALSE, params = c("pi"))


