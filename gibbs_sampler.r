######################################
## Gibbs Sampler for Spam Detection ##
######################################

# This code is a hierarchical Bayesian model that uses a Gibbs sampler to filter
# emails. The model assumes that the number of "spam" words in an email is
# Poisson distributed. So the inbox is a mixture of Poisson distributions with
# a latent variable of whether or not the email is spam. We are interested in
# the value of the latent variable for our classification.

# Reading the data
spam<-read.csv("spam.csv")
#spam<-spam[1:500,]
n<-length(spam[,1])
real<-rep(1,n)
real[spam[,1]=="spam"]<-0


# Setting parameters
set.seed(405)               # For testing reproducibility
Z0<-rbinom(n, 1, .5)        # Initial value of Z for Gibbs Sampler
beta<-1                     # Used in exponential prior
a<-1                        # Used in beta prior
b<-1                        # Used in beta Prior
m<-1000                     # Number of Gibbs Sampler iterations




# Create a list of spam words to test each email against
spam_words<-tolower(c("Free", "entry", "wkly", "comp", "win", "final", "tkts", "Text", "entry",
              "FreeMsg", "darling", "fun", "XxX", "entitled", "Update", "credit", "claim", "subscription",
              "special", "pleased", "valued", "jackpot", "prize"))

# Assigning spam word counts to each email
# Create X vector representing the number of spam words in each email
X<-rep(0,n)

# Search through each email to see update spam word count
vec<-spam_words
for(word in vec){
  for(i in 1:n){
    if(grepl(word, tolower(spam[i,2]))){
      X[i]<-X[i]+1
    }
  }
}


# Defining functions to make easy sampling
samp<-function(X, l0, l1, p){
  temp<-c()
  p.0<-(1-p)*dpois(X, l0)
  p.1<-p*dpois(X, l1)
  q<-p.1/(p.1+p.0)
  for(i in 1:n){
    temp<-c(temp, rbinom(1,1,q[i]))
  }
  return(temp)
}

# Running the Gibbs Sampler
Z<-vector(mode='list', length=m+1)                               # Creating the Z vector
Z[[1]]<-Z0                                                       # Initializing the Z vector
lambda<-vector(mode="list", length=m+1)
mix_param<-c()
for(j in 1:m){
  Y<-Z[[j]]
  l0<-rgamma(1,sum(X*(1-Y))+1, scale=beta/(beta*sum(1-Y)+1))     # Sample conditional lambda0
  l1<-rgamma(1,sum(X*Y)+1, scale=beta/(beta*sum(Y)+1))           # Sample conditional lambda1
  lambda[[j]]<-c(l0, l1)                                         # Saving iterations of the parameters
  p<-rbeta(1, a+sum(Y), b+sum(1-Y))                              # Sample conditional mixing parameter
  mix_param<-c(mix_param, p)                                     # Saving iterations of the mixing parameter
  Z[[j+1]]<-samp(X, l0, l1, p)                                   # Sample the labeling variables
  if(j%%500==0){
    print(j)
  }
}

# Finding the posterior mean of Z
Z.bar<-rep(0,n)
for(i in 1:m+1){
  Z.bar<-Z.bar+Z[[i]]
}
Z.bar<-Z.bar/m

# Rounding the estimates to binary spam/not-spam labels
Z.bar<-round(Z.bar)


# Error Analysis
# Finding the accuracy of the algorithm
error<-sum((Z.bar-real)^2)
percent_correct<-1-error/n
percent_correct                        # 91% accurate

# Counting the number of emails marked as not spam when they actually are
spam_count<-0
for(i in 1:n){
  if(Z.bar[i]==1 & real[i]==0){
    spam_count<-spam_count+1
  }
}
spam_count


# Counting the number of non-spam emails marked as spam
ham_count<-0
for(i in 1:n){
  if(Z.bar[i]==0 & real[i]==1){
    ham_count<-ham_count+1
  }
}
ham_count


# Determining the real values of the parameters of our model
l0_real<-mean(X[real==0])
l1_real<-mean(X[real==1])
p_real<-sum(real)/n

# Finding the spam word count of the mislabeled emails
missed<-X[real!=Z.bar]
l_missed<-mean(missed)                    # Average spam word count of the missed emails
l_01<-(l0_real-l1_real)/2                 # Average spam word count between the two groups
                                          # l_missed is approx. l_01, the missed emails did not give enough
                                          #    information to definitively classify them. Improving the
                                          #    spam word bank would better separate these groups


# Finding and plotting the error at each iteration of the Gibbs Sampler
error_vec<-c()
for(i in 1:(m+1)){
  error_vec<-c(error_vec, sum((Z[[i]]-real)^2))
}

plot(error_vec, xlab="Gibbs Sampler Iteration Numebr", ylab=" ", main="Error vs Iteration Number", type="l")





