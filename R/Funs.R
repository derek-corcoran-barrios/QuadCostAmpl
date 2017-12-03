#' Generates an AMPL dat file from a list of Stacks
#'
#' Generates an AMPL dat file from a Stack in which each file is the projection
#' of a species distribution model into a time slice
#' @param Stacklist a list of Stacks with the space where the species will be
#' inhabiting
#' @param Dist the maximum dispersal distance of the species modeled in the
#' stack
#' @param name the name of the .dat file that will be exported
#' @param nchains the number of chains to go through
#' @param costlayer raster with the costs of each cell
#' @return exports a .dat file to feed the AMPL model
#' @examples
#' \dontrun{
#' data("BinSpp")
#' data("Cost")
#' MultiSppQuad(Stacklist = BinSpp, Dist = 1000000, name = "Two", costlayer = Cost)
#' }
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom gdistance accCost
#' @importFrom gdistance geoCorrection
#' @importFrom gdistance transition
#' @importFrom gdistance transitionMatrix
#' @importFrom igraph E
#' @importFrom igraph "E<-"
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest.paths
#' @importFrom magrittr "%>%"
#' @importFrom Matrix cBind
#' @importFrom Matrix rBind
#' @importFrom raster ncell
#' @importFrom raster nlayers
#' @importFrom raster values
#' @importFrom raster xyFromCell
#' @importFrom tidyr unite_
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Javier Fajardo <javierfajnolla@gmail.com >
#' @export


MultiSppQuad <- function(Stacklist, Dist, name, nchains = 100, costlayer){

  Masklayer <- costlayer
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)
  for (i in 1:length(Stacklist)){
    Stacklist[[i]] <- Stacklist[[i]] * Masklayer
  }
  accCost2 <- function(x, fromCoords) {

    fromCells <- cellFromXY(x, fromCoords)
    tr <- transitionMatrix(x)
    tr <- rBind(tr, rep(0, nrow(tr)))
    tr <- cBind(tr, rep(0, nrow(tr)))
    startNode <- nrow(tr)
    adjP <- cbind(rep(startNode, times = length(fromCells)), fromCells)
    tr[adjP] <- Inf
    adjacencyGraph <- graph.adjacency(tr, mode = "directed", weighted = TRUE)
    E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
    return(shortest.paths(adjacencyGraph, v = startNode, mode = "out")[-startNode])
  }
  Suitabilities <- list()
  for(j in 1:length(Stacklist)){
    Suitability <- list()
    for (i in 1:nlayers(Stacklist[[j]])){
      temp <- data.frame(Suitability = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i-1)
      Suitability[[i]] <- temp[complete.cases(temp),]
    }
    Suitabilities[[j]]<- do.call("rbind", Suitability)
    Suitabilities[[j]]$Spp <- names(Stacklist)[j]
  }
  Suitability <- do.call("rbind", Suitabilities)
  s <- Suitability %>% group_by(ID) %>% summarise(SUMA = sum(Suitability)) %>% filter(SUMA > 0)
  Suitability <- Suitability[Suitability$ID %in% s$ID,]


  Spps <- unique(Suitability$Spp)

  Suitability <- Suitability[,c(4,1,2,3)]

  Suitabilities <- list()
  for (i in Spps){
    Suitabilities[[i]] <- dplyr::filter(Suitability, Spp == i)

    temp <-  split(Suitabilities[[i]], Suitabilities[[i]]$Time)
    Suitabilities[[i]] <- do.call(cbind, lapply(1:length(temp), function(i){
      if (i == 1){
        setNames(data.frame(paste("[",temp[[i]][["Spp"]],",",temp[[i]][["ID"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["Suitability"]),
                 c("V", paste("T", i, sep = ""), i-1))
      } else if (i == length(temp)){
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Suitability"], rep("\n", NROW(temp[[i]]))),
                 c(paste("T", i, sep = ""), i-1, "line"))
      } else {
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Suitability"]),
                 c(paste("T", i, sep = ""), i-1))
      }
    }))
  }

  Suitability <-do.call("rbind", Suitabilities)

  conns <- list()
  for(j in 1:length(Stacklist)){
    Raster <- sum(Stacklist[[j]])

    Raster[values(Raster) > 0] = 1
    Raster[values(Raster) == 0] = NA

    h16  <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)

    h16   <- geoCorrection(h16, scl=FALSE)

    ID <-c(1:ncell(Raster))[!is.na(values(Raster))]

    B <- xyFromCell(Raster, cell = ID)

    connections <- list()
    #For each pair of cells in B
    for (i in 1:nrow(B)){
      #Create a temporal raster for each row with the distance from cell xy to all other cells
      temp <- accCost2(h16,B[i,])
      index <- which(temp < Dist)
      connections[[i]] <- cbind(ID[i], index, temp[index])
    }
    #Get everything together as a large data frame
    connections <- do.call("rbind", connections)
    connections <- as.data.frame(connections)
    colnames(connections) <- c("from", "to", "dist")
    connections$Sp <- names(Stacklist)[j]
    conns[[j]] <- connections
  }
  connections <- conns
  connections <- do.call("rbind", connections)

  Nchains <- data.frame(Spp = Spps, Nchains = nchains, Space = "\n")
  Cost <- values(costlayer)[unique(unique(connections$to), unique(connections$to))]


  sink(paste0(name, ".dat"))
  cat(c("set V :=", unique(unique(connections$to), unique(connections$to)), ";"))
  cat("\n")
  cat("\n")
  cat(c("set SP :=", names(Stacklist), ";"))
  cat("\n")
  cat("\n")
  cat(c("set c :=", Cost, ";"))
  cat("\n")
  cat("\n")
  cat(c("set E :=", paste0("(",unique(unite_(connections, col = "V", sep = ",", from = c("from", "to"))$V), ")"), ";"))
  cat("\n")
  cat("\n")
  cat(paste0("param T:= ", (nlayers(Stacklist[[1]])-1),";"))
  cat("\n")
  cat("\n")
  cat("param u :=")
  cat("\n")
  cat(do.call(paste, Suitability))
  cat(";")
  cat("\n")
  cat("param nchains := ")
  cat("\n")
  cat(do.call(paste, Nchains))
  cat(";")
  cat("\n")
  sink()
  return(list(connections = connections, Suitability = Suitability))
}
