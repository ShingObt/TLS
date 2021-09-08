# .pol create polygons as regular squares overlapped 0.05

.pol <- function(data, dist){
  
  .xpol <- c(data$x.1-dist, data$x.2+dist, data$x.3+dist, data$x.4-dist)
  .ypol <- c(data$y.1-dist, data$y.2-dist, data$y.3+dist, data$y.4+dist)
  
  return(sp::Polygons(list(sp::Polygon(cbind(.xpol,.ypol))), ID = data$id))
  
}
.ncr.remove.slice.double <- function(data){
  
  # Select necessary fields from original txt file of point cloud
  
  .data <- data[,c("point", "x", "y", "z")]
  
  
  # Create x and y coordinates for grid
  
  .x = seq(min(.data$x), max(.data$x))
  .y = seq(min(.data$y), max(.data$y))
  
  # Empty data frame where coordinates neccesaries for
  # creating grid will be saved
  
  .grid <- data.frame(id = as.character(),
                      x.1 = as.numeric(), x.2 = as.numeric(),
                      x.3 = as.numeric(), x.4 = as.numeric(),
                      y.1 = as.numeric(), y.2 = as.numeric(),
                      y.3 = as.numeric(), y.4 = as.numeric())
  
  
  # Fill the empty data frame .grid created just before
  
  for (i in 1:(length(.x)-1)) {
    for (j in 1:(length(.y)-1)) {
      
      .out <- data.frame(id = paste("x", i, "y", j, sep="."),
                         x.1 = .x[i], x.2 = .x[i+1], x.3 = .x[i+1], x.4 = .x[i],
                         y.1 = .y[j], y.2 = .y[j], y.3 = .y[j+1], y.4 = .y[j+1])
      
      .grid <- rbind(.grid, .out)
      
    }
  }
  
  .row.names <- .grid$id
  
  # Split .grid by id into a list for using lapply function and
  # applying .pol function over all elements of the list
  
  .grid <- split(.grid, .grid$id)
  
  # Apply .pol function to every element of the list and create
  # and object SpatialPolygons
  
  .grid <- sp::SpatialPolygons(lapply(.grid, .pol, dist = 0.05))
  
  
  # Generate and SpatialPoints object to extract those points
  # which are avor the polygons of .grid
  
  .pts <- .data[, c("x", "y")]
  # dimnames(.pts)[[1]] <- c(.data$point)
  .pts <- sp::SpatialPoints(.pts)
  
  .attributes <- .data[, c("point", "x", "y", "z")]
  
  .pts = sp::SpatialPointsDataFrame(.pts, .attributes)
  
  .attributes <- data.frame(row.names = .row.names)
  
  .grid = sp::SpatialPolygonsDataFrame(.grid, .attributes)
  
  .pts <- sp::over(.grid, .pts, returnList = TRUE)
  
  .dat <- lapply(.pts, as.matrix, ncol = 4)
  .dat <- .dat[names(which(lapply(.dat, length) > 4))]
  # .dat <- .dat[names(which(lapply(.dat, length) < 40000))]
  
  
  .ncr <- do.call(rbind, lapply(.dat, ncr_point_cloud_double))
  
  .ncr <- .ncr[which(.ncr$ncr > 0 & .ncr$ncr < 9999), ]
  
  .data <- merge(data, .ncr, by = "point", all = TRUE)
  
  .data <- .data[!duplicated(.data), ]
  
  return(.data)
  
}

















treeDetection=function (data, dbh.min = 7.5, dbh.max = 200, ncr.threshold = 0.1, 
                        tls.resolution = list(point.dist = 7.67, tls.dist = 10), breaks = c(1, 1.3, 1.6), plot.attributes = NULL, 
          save.result = TRUE, dir.result = NULL) 
{
  if (is.null(dir.result)) 
    dir.result <- getwd()
  .dbh.min <- dbh.min/100
  .dbh.max <- dbh.max/100
  .point.dist <- tls.resolution$point.dist/1000
  .tls.dist <- tls.resolution$tls.dist
  .vertical.angle <- tls.resolution$vertical.angle * 2 * pi/360
  .horizontal.angle <- tls.resolution$horizontal.angle * 2 * 
    pi/360
  print(tls.resolution)
  if (is.null(.point.dist)) {
    .alpha.v <- .vertical.angle
    .alpha.h <- .horizontal.angle
  }else{
    .alpha.v <- atan((.point.dist/2)/.tls.dist) * 2
    .alpha.h <- .alpha.v
  }
  .filteraux <- data.frame(cluster = as.numeric(), center.x = as.numeric(), 
                           center.y = as.numeric(), center.phi = as.numeric(), center.rho = as.numeric(), 
                           center.r = as.numeric(), center.theta = as.numeric(), 
                           radius = as.numeric(), num.points = as.numeric(), num.points.hom = as.numeric(), 
                           phi.left = as.numeric(), phi.right = as.numeric(), arc.circ = as.numeric(), 
                           sec = as.numeric())
  for (cuts in breaks) {
    message("Computing section: ", cuts, " m")
    .cut <- data[which(data$z > (cuts - 0.1) & data$z < (cuts + 
                                                           0.1)), , drop = FALSE]
    .cut <- .ncr.remove.slice.double(.cut)
    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), 
                 , drop = FALSE]
    .cut <- .cut[which(.cut$z > (cuts - 0.05) & .cut$z < 
                         (cuts + 0.05)), , drop = FALSE]
    .eps <- .dbh.min/2
    .dbscan <- dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], 
                              eps = .eps)
    .cut$cluster <- .dbscan$cluster
    .cut <- .cut[which(.cut$cluster > 0), , drop = FALSE]
    if (nrow(.cut) < 1) 
      next
    .cut$sec <- cuts
    .pb <- progress::progress_bar$new(total = length(unique(.cut$cluster)))
    .filter <- data.frame(cluster = as.numeric(), center.x = as.numeric(), 
                          center.y = as.numeric(), center.phi = as.numeric(), 
                          center.rho = as.numeric(), center.r = as.numeric(), 
                          center.theta = as.numeric(), radius = as.numeric(), 
                          num.points = as.numeric(), num.points.hom = as.numeric(), 
                          phi.left = as.numeric(), phi.right = as.numeric(), 
                          arc.circ = as.numeric(), occlusion = as.numeric())
    for (.i in unique(.cut$cluster)) {
      .pb$tick()
      .dat <- .cut[which(.cut$cluster == .i), , drop = FALSE]
      .x.rang <- max(.dat$x) - min(.dat$x)
      .y.rang <- max(.dat$y) - min(.dat$y)
      .phi.rang <- max(.dat$phi) - min(.dat$phi)
      .rho.rang <- max(.dat$rho) - min(.dat$rho)
      .x.cent <- (.x.rang/2) + min(.dat$x)
      .y.cent <- (.y.rang/2) + min(.dat$y)
      .phi.cent <- (.phi.rang/2) + min(.dat$phi)
      .rho.cent <- (.rho.rang/2) + min(.dat$rho)
      .ancho.malla <- (max(.x.rang, .y.rang)/2) * 1.5
      .xmin <- .x.cent - .ancho.malla
      .ymin <- .y.cent - .ancho.malla
      .xmax <- .x.cent + .ancho.malla
      .ymax <- .y.cent + .ancho.malla
      .ancho.malla.2 <- (max(.phi.rang, .rho.rang)/2)
      .phimin <- .phi.cent - .ancho.malla.2
      .rhomin <- .rho.cent - .ancho.malla.2
      .phimax <- .phi.cent + .ancho.malla.2
      .rhomax <- .rho.cent + .ancho.malla.2
      .x.values <- seq(from = .xmin, to = .xmax, by = 0.03)
      .y.values <- seq(from = .ymin, to = .ymax, by = 0.03)
      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
      for (.i in 1:length(.x.values)) {
        for (.j in 1:length(.y.values)) {
          .den <- .dat[which(.dat$x <= ((.x.values[.i]) + 
                                          0.015) & .dat$x > ((.x.values[.i]) - 0.015) & 
                               .dat$y <= ((.y.values[.j]) + 0.015) & .dat$y > 
                               ((.y.values[.j]) - 0.015)), , drop = FALSE]
          .density[.j, .i] <- ifelse(nrow(.den) < 1, 
                                     NA, nrow(.den))
        }
      }
      .threeshold <- mean(.density, na.rm = T)
      if (is.nan(.threeshold)) {
        next
      }
      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
      .remove <- data.frame(point = as.numeric())
      for (.i in 1:length(.x.values)) {
        for (.j in 1:length(.y.values)) {
          .den <- .dat[which(.dat$x <= ((.x.values[.i]) + 
                                          0.015) & .dat$x > ((.x.values[.i]) - 0.015) & 
                               .dat$y <= ((.y.values[.j]) + 0.015) & .dat$y > 
                               ((.y.values[.j]) - 0.015)), , drop = FALSE]
          .density[.j, .i] <- ifelse(nrow(.den) < 1, 
                                     NA, nrow(.den))
          if (nrow(.den) > .threeshold) {
            .rem <- data.frame(point = .den$point)
            .remove <- rbind(.remove, .rem)
          }
        }
      }
      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)
      if (nrow(.dat) < 1) {
        next
      }
      if (is.nan(mean(.dat$slope, na.rm = TRUE))) {
        .n <- (0.1/(tan(.alpha.v/2) * (mean(.dat$r)/cos(mean(.cut$slope, 
                                                             na.rm = TRUE))) * 2))
      }
      else {
        .n <- (0.1/(tan(.alpha.v/2) * (mean(.dat$r)/cos(mean(.dat$slope, 
                                                             na.rm = TRUE))) * 2))
      }
      if (mean(.cut$slope, na.rm = TRUE) > 0.5) {
        .n <- 0.7 * .n
      }
      .x2.values <- seq(from = min(.dat$phi), to = max(.dat$phi), 
                        by = .alpha.h)
      .density <- vector(length = length(.x2.values))
      .remove <- data.frame(point = as.numeric())
      for (.i in 1:length(.x2.values)) {
        .den <- .dat[which(.dat$phi <= ((.x2.values[.i]) + 
                                          (.alpha.h/2)) & .dat$phi > ((.x2.values[.i]) - 
                                                                        (.alpha.h/2))), ]
        .density[.i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))
        if (nrow(.den) > 1) {
          .rem <- data.frame(point = .den$point)
          .remove <- rbind(.remove, .rem)
        }
      }
      .density <- ifelse(is.nan(.density), NA, .density)
      if (is.nan(mean(.density, na.rm = TRUE))) {
        next
      }
      if (max(.density[which(!is.na(.density))], na.rm = T) < 
          floor(.n)) {
        next
      }
      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)
      if (nrow(.dat) < 1) {
        next
      }
      .num.points <- nrow(.dat)
      .num.points.hom <- nrow(.dat[which(.dat$prob.selec == 
                                           1), , drop = FALSE])
      .x.values <- seq(from = .xmin, to = .xmax, by = 0.01)
      .y.values <- seq(from = .ymin, to = .ymax, by = 0.01)
      .matriz <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
      for (.i in 1:length(.x.values)) {
        for (.j in 1:length(.y.values)) {
          .variance <- stats::var(raster::pointDistance(cbind(.dat$x, 
                                                              .dat$y), c(.x.values[.i], .y.values[.j]), 
                                                        lonlat = FALSE))
          .matriz[.j, .i] <- .variance
        }
      }
      .a <- which(.matriz == min(.matriz), arr.ind = TRUE)
      .center.x <- .x.values[.a[2]]
      .center.y <- .y.values[.a[1]]
      .center.phi <- atan2(.center.y, .center.x)
      .center.phi <- ifelse(.center.phi < 0, .center.phi + 
                              (2 * pi), .center.phi)
      .center.rho <- sqrt(.center.x^2 + .center.y^2)
      .center.r <- sqrt(.dat$sec[1]^2 + .center.rho^2)
      .center.theta <- atan2(.dat$sec[1], .center.rho)
      .radio <- mean(raster::pointDistance(cbind(.dat$x, 
                                                 .dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), 
                                           lonlat = FALSE))
      .dat$dist <- raster::pointDistance(cbind(.dat$x, 
                                               .dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), 
                                         lonlat = FALSE)
      if (stats::quantile(.dat$rho, prob = 0.05) > .center.r) {
        next
      }
      .dat <- .dat[order(.dat$dist, decreasing = FALSE), 
                   , drop = FALSE]
      if (stats::quantile(.dat$dist, prob = 0.05) < (.radio/2)) {
        next
      }
      if ((max(.dat$phi) - min(.dat$phi)) < pi) {
        .dat.2 <- .dat[order(.dat$phi, decreasing = F), 
                       , drop = FALSE]
      }
      else {
        .dat.2 <- .dat
        .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + 
                               (2 * pi), .dat.2$phi)
        .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), 
                         , drop = FALSE]
      }
      .ratio <- nrow(.dat.2)/(.n * ((max(.dat.2$phi) - 
                                       min(.dat.2$phi))/.alpha.h))
      if (.ratio < 0.5) {
        next
      }
      .pto.left <- stats::quantile(.dat.2$phi, prob = 0.01)
      .rho.left <- mean(.dat.2$rho[which(.dat.2$phi <= 
                                           .pto.left)])
      .phi.left <- mean(.dat.2$phi[which(.dat.2$phi <= 
                                           .pto.left)])
      .pto.right <- stats::quantile(.dat.2$phi, prob = 0.99)
      .rho.right <- mean(.dat.2$rho[which(.dat.2$phi >= 
                                            .pto.right)])
      .phi.right <- mean(.dat.2$phi[which(.dat.2$phi >= 
                                            .pto.right)])
      .phi.cent <- max(.dat.2$phi) - ((max(.dat.2$phi) - 
                                         min(.dat.2$phi))/2)
      .rho.cent <- mean(.dat.2$rho[which(round(.dat.2$phi, 
                                               3) >= round(.phi.cent - .alpha.h, 3) & round(.dat.2$phi, 
                                                                                            3) <= round(.phi.cent + .alpha.h, 3))])
      if (is.nan(.rho.cent)) {
        next
      }
      .arc.circ <- ifelse(.rho.left > .rho.cent & .rho.right > 
                            .rho.cent, 1, 0)
      .phi.left <- ifelse(.phi.left > (2 * pi), .phi.left - 
                            (2 * pi), .phi.left)
      .phi.right <- ifelse(.phi.right > (2 * pi), .phi.right - 
                             (2 * pi), .phi.right)
      .dat.2$n <- c(1:nrow(.dat.2))
      .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$phi, 
                                  method = "pearson"), silent = TRUE)
      if (class(.cor) == "try-error") {
        next
      }
      else {
        .occlusion <- .cor[[4]]
      }
      .cv <- stats::sd(raster::pointDistance(cbind(.dat$x, 
                                                   .dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), 
                                             lonlat = FALSE))/.radio
      if (.cv > 0.1) {
        next
      }
      .n.w.ratio <- stats::sd(.dat$z)/sqrt(stats::sd(.dat$x)^2 + 
                                             stats::sd(.dat$y)^2)
      if (.n.w.ratio > 1.5) {
        next
      }
      .salida <- data.frame(cluster = .dat$cluster[1], 
                            center.x = .center.x, center.y = .center.y, center.phi = .center.phi, 
                            center.rho = .center.rho, center.r = .center.r, 
                            center.theta = .center.theta, radius = .radio, 
                            num.points = .num.points, num.points.hom = .num.points.hom, 
                            phi.left = .phi.left, phi.right = .phi.right, 
                            arc.circ = .arc.circ, occlusion = .occlusion)
      .filter <- rbind(.filter, .salida)
    }
    .filter$tree <- ifelse(.filter$arc.circ == 1, 1, ifelse(.filter$arc.circ == 
                                                              0 & .filter$occlusion > 0.995, 1, 0))
    .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]
    .filter$tree <- ifelse(.filter$radius > (.dbh.max/2) | 
                             .filter$radius < (.dbh.min/2), 0, 1)
    .filter <- subset(.filter, .filter$tree == 1)
    if (nrow(.filter) < 1) {
      .filter1.0 <- data.frame(cluster = as.numeric(), 
                               center.x = as.numeric(), center.y = as.numeric(), 
                               center.phi = as.numeric(), center.rho = as.numeric(), 
                               center.r = as.numeric(), center.theta = as.numeric(), 
                               radius = as.numeric(), num.points = as.numeric(), 
                               num.points.hom = as.numeric(), phi.left = as.numeric(), 
                               phi.right = as.numeric(), arc.cir = as.numeric(), 
                               sec = as.numeric())
    }
    else {
      .filter1.0 <- .filter[, c("cluster", "center.x", 
                                "center.y", "center.phi", "center.rho", "center.r", 
                                "center.theta", "radius", "num.points", "num.points.hom", 
                                "phi.left", "phi.right", "arc.circ"), drop = FALSE]
      .filter1.0$sec <- cuts
    }
    .filteraux <- rbind(.filteraux, .filter1.0)
  }
  .filter <- .filteraux
  if (nrow(.filter) < 1) 
    stop("No tree was detected")
  .dbscan <- dbscan::dbscan(.filter[, c("center.x", "center.y"), 
                                    drop = FALSE], eps = max(.filter$radius), minPts = 1)
  .filter$cluster <- .dbscan$cluster
  .filter <- .filter[order(.filter$cluster, .filter$sec), , 
                     drop = FALSE]
  .taper <- .filter[, c("cluster", "sec", "radius"), drop = FALSE]
  .lm <- stats::lm(radius ~ sec, data = .taper)
  .slope <- stats::coef(.lm)[2]
  .filter$dif <- 1.3 - .filter$sec
  .filter$radio.est <- ifelse(.filter$dif == 0, .filter$radius, 
                              .filter$radius + .slope * .filter$dif)
  .filter$radio.est <- ifelse(is.na(.filter$radio.est), .filter$radius, 
                              .filter$radio.est)
  .radio.est <- tapply(.filter$radio.est, .filter$cluster, 
                       mean)
  .tree <- data.frame(tree = tapply(.filter$cluster, .filter$cluster, 
                                    mean, na.rm = TRUE), center.x = tapply(.filter$center.x, 
                                                                           .filter$cluster, mean, na.rm = TRUE), center.y = tapply(.filter$center.y, 
                                                                                                                                   .filter$cluster, mean, na.rm = TRUE), center.phi = tapply(.filter$center.phi, 
                                                                                                                                                                                             .filter$cluster, mean, na.rm = TRUE), center.rho = tapply(.filter$center.rho, 
                                                                                                                                                                                                                                                       .filter$cluster, mean, na.rm = TRUE), center.r = tapply(.filter$center.r, 
                                                                                                                                                                                                                                                                                                               .filter$cluster, mean, na.rm = TRUE), center.theta = tapply(.filter$center.theta, 
                                                                                                                                                                                                                                                                                                                                                                           .filter$cluster, mean, na.rm = TRUE), horizontal.distance = tapply(.filter$center.rho, 
                                                                                                                                                                                                                                                                                                                                                                                                                                              .filter$cluster, mean, na.rm = TRUE), radius = .radio.est, 
                      phi.left = tapply(.filter$phi.left, .filter$cluster, 
                                        mean, na.rm = TRUE), phi.right = tapply(.filter$phi.right, 
                                                                                .filter$cluster, mean, na.rm = TRUE), partial.occlusion = tapply(.filter$arc.circ, 
                                                                                                                                                 .filter$cluster, mean, na.rm = TRUE), num.points = tapply(.filter$num.points, 
                                                                                                                                                                                                           .filter$cluster, mean, na.rm = TRUE), num.points.hom = tapply(.filter$num.points.hom, 
                                                                                                                                                                                                                                                                         .filter$cluster, mean, na.rm = TRUE))
  .tree$partial.occlusion <- ifelse(.tree$partial.occlusion == 
                                      0, 1, 0)
  .tree$dbh <- .tree$radius * 200
  .filter$filter <- ifelse(.filter$sec == 1.3 & .filter$arc.circ == 
                             1, 1, 0)
  .filter2 <- subset(.filter, .filter$filter == 1)
  if (nrow(.filter2) < 1) 
    .filter2 <- .filter
  .filter2$points.radio <- .filter2$num.points/.filter2$radio
  .filter2$points.radio.hom <- .filter2$num.points.hom/.filter2$radio
  .tree$points.m <- mean(.filter2$points.radio)
  .tree$points.m.hom <- mean(.filter2$points.radio.hom)
  .tree$num.points.est <- .tree$points.m * .tree$radius
  .tree$num.points.hom.est <- .tree$points.m.hom * .tree$radius
  if (is.null(data$id)) {
    .tree <- .tree[, c("tree", "center.x", "center.y", "center.phi", 
                       "phi.left", "phi.right", "horizontal.distance", "dbh", 
                       "num.points", "num.points.hom", "num.points.est", 
                       "num.points.hom.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "phi.left", 
                         "phi.right", "horizontal.distance", "dbh", "num.points", 
                         "num.points.hom", "num.points.est", "num.points.hom.est", 
                         "partial.occlusion")
  }
  else {
    .tree$id <- data$id[1]
    .tree$file <- data$file[1]
    .tree <- .tree[, c("id", "file", "tree", "center.x", 
                       "center.y", "center.phi", "phi.left", "phi.right", 
                       "horizontal.distance", "dbh", "num.points", "num.points.hom", 
                       "num.points.est", "num.points.hom.est", "partial.occlusion"), 
                   drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", 
                         "phi", "phi.left", "phi.right", "horizontal.distance", 
                         "dbh", "num.points", "num.points.hom", "num.points.est", 
                         "num.points.hom.est", "partial.occlusion")
  }
  if (!is.null(plot.attributes)) 
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)
  if (isTRUE(save.result)) {
    utils::write.csv(.tree, file = file.path(dir.result, 
                                             "tree.list.tls.csv"), row.names = FALSE)
  }
  return(.tree)
}




















TLSnormalize <- function(las,
                      max.dist = NULL, min.height = NULL, max.height = NULL,
                      algorithm.dtm = "tin", res.dtm = 0.2,
                      id = NULL, file=NULL,
                      dir.data = NULL, save.result = TRUE, dir.result = NULL){
  
  .pb <- progress::progress_bar$new(total = 11)
  .pb$tick()
  
  
  # Obtaining working directory for loading files
  if(is.null(dir.data))
    dir.data <- getwd()
  
  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()
  
  
  # Loading input (LAS file)
  
  .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyz")))
  
  .pb$tick()
  
  
  # Giving the same scale factor to all coordinates
  
  .las@header@PHB[["X scale factor"]] <- 0.001
  .las@header@PHB[["Y scale factor"]] <- 0.001
  .las@header@PHB[["Z scale factor"]] <- 0.001
  
  
  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane
  
  if(!is.null(max.dist)) {
    
    .las <- lidR::clip_circle(.las, 0, 0, max.dist)
    
  }
  
  .pb$tick()
  
  
  # Normalize
  
  # .ws  <- seq(3,12, 4)
  # .th  <- seq(0.1, 1.5, length.out = length(.ws))
  # .data <- lidR::classify_ground(.las, algorithm = lidR::pmf(.ws, .th), last_returns = FALSE)
  
  .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(), last_returns = FALSE)))
  
  .pb$tick()
  
  
  # Generaion of Digital Terrain Model (DTM)
  
  
  if(algorithm.dtm == "knnidw")
    .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::knnidw())))
  
  
  if(algorithm.dtm == "tin")
    .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::tin())))
  
  
  .dtm[.dtm < min(.data@data$Z)] <- NA
  
  .pb$tick()
  
  # Estimating slope
  
  raster::crs(.dtm) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83"
  .slope <- raster::terrain(.dtm, opt=c('slope'), unit='radians')
  
  .pb$tick()
  
  if(mean(.slope@data@values, na.rm = TRUE) > 0.3){
    
    
    # Normalize
    
    .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(sloop_smooth = TRUE), last_returns = FALSE)))
    
    
    # Generaion of Digital Terrain Model (DTM)
    
    if(algorithm.dtm == "knnidw")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::knnidw())))
    
    
    if(algorithm.dtm == "tin")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::tin())))
    
    
    .dtm[.dtm < min(.data@data$Z)] <- NA
    
    
  }
  
  .pb$tick()
  
  
  # Normalization of cloud data
  
  .data <- suppressWarnings(suppressMessages(lidR::normalize_height(.data, .dtm, add_lasattribute = FALSE, na.rm = TRUE)))
  
  .pb$tick()
  
  # Assigning slope to point cloud
  
  .data <- lidR::merge_spatial(.data, .slope, "slope")
  
  .pb$tick()
  
  # Test for getting a smaller file - data frame
  
  .data <- data.frame(.data@data)
  
  
  # Removing points classified as ground
  
  .data <- subset(.data, .data$Classification == 1)
  
  
  # Extracting coordinates values
  
  .data <- .data[, c("X", "Y", "Z", "slope"), drop = FALSE]
  colnames(.data) <- c("x", "y", "z", "slope")
  
  
  # Low point filtering
  
  if(!is.null(min.height)) {
    
    .data <- .data[which(.data$z > min.height), , drop = FALSE]
    
  }
  
  
  # High point filtering
  
  if(!is.null(max.height)) {
    
    .data <- .data[which(.data$z < max.height), , drop = FALSE]
    
  }
  
  # Transformation to other coordinate systems
  
  # Cylindrical coordinate system (https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
  # rho, axial distance or radial distance (euclidean distance from the z-axis to the point P)
  # phi, azimuth is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane
  # z, axial coordinate or height z is the signed distance from the chosen plane to the point P
  .data$rho <- sqrt(.data$x ^ 2 + .data$y ^ 2)
  .data$phi <- atan2(.data$y, .data$x)
  .data$phi <- ifelse(.data$phi < 0, .data$phi + (2 * pi), .data$phi)
  
  # Spherical coordinates system (https://en.wikipedia.org/wiki/Spherical_coordinate_system)
  # r, radius or radial distance is the Euclidean distance from the origin O to P
  # theta, inclination (or polar angle) is the angle between the zenith direction and the line segment OP
  # phi, azimuth is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane
  
  .data$r <- sqrt(.data$z ^ 2 + .data$rho ^ 2)
  .data$theta <- atan2(.data$z, .data$rho)
  
  .data$point <- (1:nrow(.data))
  
  
  # Point crooping process
  # This is a previous step to obtain a homogeneous density of points in the space
  # This is based on the principle that closer objects (with the same size and shape)
  # have more probability to recieve points
  
  .data$prob <- (.data$r / max(.data$r)) ^ 2
  .data$prob.random <- stats::runif(nrow(.data))
  .data$prob.selec <- ifelse(.data$prob >= .data$prob.random, 1, 0)
  
  
  # Assign id
  
  if(!is.null(id)){
    
    .data$id <- id
    
  } else {
    
    .data$id <- 1
    
  }
  
  
  # File name
  
  if(!is.null(file)){
    
    .data$file <- file
    
  } else {
    
    .data$file <- paste(.data$id[1], ".txt", sep = "")
    
  }
  
  
  .data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "prob", "prob.selec"), drop = FALSE]
  
  .pb$tick()
  
  # Saving data
  
  # Obtaining working directory for saving files
  
  if(isTRUE(save.result)){
    
    .data.red <- .data[which(.data$prob.selec == 1), , drop = FALSE]
    
    vroom::vroom_write(.data.red, file = file.path(dir.result, .data.red$file[1]), delim = ",", progress = FALSE)
    
  }
  
  # .pb$tick()
  
  return(.data)
  
  
  
}
