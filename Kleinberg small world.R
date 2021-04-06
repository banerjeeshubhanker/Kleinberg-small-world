#Loading the library
library(igraph)
#The size of one dimension in the 
#2D garph 
length<-7


#Make two directed lattices, one for Watts small world
#and one for Kleinberg's small world
lat_watts <- make_lattice( c(length,length), directed = TRUE,mutual = TRUE, circular = T )
lat_kleinberg<-lat_watts
index<-0

#Assign attributes to each vertex in the graphs
for (i in 1:length) {
  for (j in 1:length) {
    index<-index+1
    V(lat_watts)[index]$name<-index
    V(lat_watts)[index]$locname<-paste0(i,",",j)
    V(lat_watts)[index]$location<-list(c(i,j))
  }
}
index<-0
for (i in 1:length) {
  for (j in 1:length) {
    index<-index+1
    V(lat_kleinberg)[index]$name<-index
    V(lat_kleinberg)[index]$locname<-paste0(i,",",j)
    V(lat_kleinberg)[index]$location<-list(c(i,j))
  }
}

#Value of rewiring parameter
beta=0.08


#Distributed Search algorithm
millgram_search <- function(g, s, t, path){
  
  #Getting the source and target vertices
  source<-V(g)[s]
  target<-V(g)[t]
  
  #Checks if it has reached the target/destination
  if(source$name==target$name){
    return (path)
  }
  
  #Finding the neighbors of the current node
  neighbours <- neighbors(g, source, mode = "out")
  #Unvisited nodes of the current node
  neighbours<- neighbours [! neighbours %in% path]
  
  #If no more neighbors left then return
  if(length(neighbours)==0){
    return (path)
  }
  
  #Intialize distances to infinity
  dist=Inf
  closest_neighb = NULL
  
  #If any of the neighbors is a target node then return path
  for(i in 1:length(neighbours)){
    if(neighbours[i]$name==target$name){
      path<-c(path,target$name)
      return (path)
    }
    
    #Finding the closest neighbor on the basis of manhattan distance
    dist1<-manhattan(target,neighbours[i])
    if (dist1<dist){
      dist<-dist1
      closest_neighbour<-neighbours[i]
    }
  }
  
  #Greedy approach, go to the neighbor closest to the target
  path<-c(path,closest_neighbour$name)
  return (millgram_search(g,closest_neighbour$name,target, path))
}

#Manhattan distance between two vertices in the graph
manhattan<- function(v1,v2){
  a<-unlist(v1$location)
  b<-unlist(v2$location)
  
  dx <- min(abs(a[1]-b[1]),length-abs(a[1]-b[1]))
  dy <- min(abs(a[2]-b[2]), length-abs(a[2]-b[2]));
  
  
  return(dx+dy)
}

#There are 49 nodes in the current lattice therefor 49*48 paths which 
#can be approximated as 49*49=2401 distances
Kleinberg_matrix<-matrix(1:2401, nrow = 49, ncol = 49)


#Calculating the ln N/distance(u,w)^2, the
#kleinberg factors, K = ln 49
normalizing_factor<-log(49)
for(i in c(1:49)){
  for(j in c(1:49)){
    if(i==j){
      Kleinberg_matrix[i, j] <- 0
    }
    else{
      distance<-manhattan(V(lat_kleinberg)[i],V(lat_kleinberg)[j])
      Kleinberg_matrix[i,j]<-normalizing_factor*(1/(distance**2))
    }
  }
}


vertices_vec<-c(1:49)

#Iterate over each of the vertices
for(i in vertices_vec){
  
  #Find all neighbors of current vertex
  neighbours1 <- neighbors(lat_kleinberg, i, mode = "out")
  
  #Filtering out all vertices of graph other than current node
  vertices_vec1<-vertices_vec[ vertices_vec != i ]
  
  #The Kleinberg factors of the current node
  Kleinberg_factors<-Kleinberg_matrix[i,vertices_vec1]
  
  #Start iterating over the neighbors
  for(j in neighbours1){
    
    #sample a vertex from the remaining vertices in vertices_vec1. This sampling is 
    #according to the distribution proposed by Kleinberg
    f<-sample(vertices_vec1, 1, prob = Kleinberg_factors,replace=F)
    
    #If edge already exists then don't rewire and go to the next neighbour
    if(lat_kleinberg[i,f] == 1 || lat_kleinberg[f,i]==1){
      next
    }
    
    #Remove edge with the current neighbor
    lat_kleinberg[i,j]<-0
    lat_kleinberg[j,i]<-0
    
    #If removing this edge and adding an edge to the sampled node 
    #creates a loop then don't delete edge with current neighbor and don't rewire
    if(length(shortest_paths(lat_kleinberg, from =
                             V(lat_kleinberg)[i], to = V(lat_kleinberg)[f], weights = NA, mode="all", 
                             output = "both")$vpath)>1){
      lat_kleinberg[i,j]<-1
      lat_kleinberg[j,i]<-1
    }
    #Else create new edge
    else{
      lat_kleinberg[i,f]<-1
      lat_kleinberg[f,i]<-1
    }
  }
}

#Perform millgram distributed search on the graph
max_distance<-c()
N <-vcount(lat_kleinberg)
distances <- c()
for(i in 1: N){
  for(j in N : 1){
    if(i!=j){
      path1<-millgram_search(lat_kleinberg,i,j, path=c(i))
      if(path1[length(path1)]==j){
        distances<-c(distances,length(path1)-1)
        if(length(path1)>length(max_distance)){
          max_distance<-path1
        }
      }
      path2<-millgram_search(lat_kleinberg,j,i, path=c(j))
      if(path2[length(path2)]==i){
        distances<-c(distances,length(path2)-1)
        if(length(path2)>length(max_distance)){
          max_distance<-path2
        }
      }
    }
  }
}  
# Aaverage path length for Kleinberg with q = 2 from distributed search
apl_kleinberg_ds<- mean(distances)
#Average path length on the basis of Graph allgorithm Djikstra to normalize
apl_kleinberg_sp<-mean_distance(lat_kleinberg,directed = T)
#Normalized Result
print(apl_kleinberg_ds/apl_kleinberg_sp)



#Calculating the Beta/distance(u,w)^0, the
#kleinberg factors for Watts-Strogatz model, Beta = 0.08
for(i in c(1:49)){
  for(j in c(1:49)){
    if(i==j){
      Kleinberg_matrix[i, j] <- 0
    }
    else{
      distance<-manhattan(V(lat_watts)[i],V(lat_watts)[j])
      Kleinberg_matrix[i,j]<-beta*(1/(distance**0))
    }
  }
}

#There are 49 nodes in the current lattice therefor 49*48 paths which 
#can be approximated as 49*49=2401 distances
vertices_vec<-c(1:49)
for(i in vertices_vec){
  neighbours1 <- neighbors(lat_watts, i, mode = "out")
  vertices_vec1<-vertices_vec[ vertices_vec != i ]
  Kleinberg_factors<-Kleinberg_matrix[i,vertices_vec1]
  for(j in neighbours1){
    f<-sample(vertices_vec1, 1, prob = Kleinberg_factors,replace=F)
    if(lat_watts[i,f] == 1 || lat_watts[f,i]==1){
      next
    }
    lat_watts[i,j]<-0
    lat_watts[j,i]<-0
    if(length(shortest_paths(lat_watts, from = V(lat_watts)[i], to = V(lat_watts)[f], weights = NA, mode="all", output = "both")$vpath)>1){
      lat_watts[i,j]<-1
      lat_watts[j,i]<-1
    }
    else{
      lat_watts[i,f]<-1
      lat_watts[f,i]<-1
    }
  }
}



max_distance<-c()

N <-vcount(lat_watts)


distances <- c()
for(i in 1: N){
  for(j in N : 1){
    if(i!=j){
      path1<-millgram_search(lat_watts,i,j, path=c(i))
      if(path1[length(path1)]==j){
        distances<-c(distances,length(path1)-1)
        if(length(path1)>length(max_distance)){
          max_distance<-path1
        }
      }
      path2<-millgram_search(lat_watts,j,i, path=c(j))
      if(path2[length(path2)]==i){
        distances<-c(distances,length(path2)-1)
        if(length(path2)>length(max_distance)){
          max_distance<-path2
        }
      }
    }
  }
} 

#Distributed search average path length 
apl_watts_ds<- mean(distances)

#Average path length using Djikstra
apl_watts_sp<-mean_distance(lat_watts,directed = T)

#Normalization
print(apl_watts_ds/apl_watts_sp)
