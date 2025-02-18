michaelisRaw <- function(mm, filname = NULL) {
  model.drm <- drm(v ~ S, data = mm, fct = MM.2())
  mml <- data.frame(S = seq(0, max(mm$S), length.out = 100))
  mml$v <- predict(model.drm, newdata = mml)
  
  bestfit <- nls(michaelis, mm, start=list(Vmax=6,Km=15))
  vmax <- summary(bestfit)$coefficients[1,1]
  Km <- summary(bestfit)$coefficients[2,1]
  
  
  vmax.std <- summary(bestfit)$coefficients[1,2]
  Km.std <- summary(bestfit)$coefficients[2,2]
  
  # Plot michaelis
  pl <- ggplot(mm, aes(x = S, y = v)) +
    theme_bw() +
    xlab("Substrate concentration") +
    ylab("Rate") +
    ggtitle(paste0("Michaelis-Menten kinetics of ", mm$enzymes[1])) +
    geom_point(alpha = 0.5) +
    geom_line(data = mml, aes(x = S, y = v), colour = "red") +
    geom_text(x = 100, y = 150, label = paste0("Vmax =", round(vmax, 2), " ± ", round(vmax.std, 2))) +
    geom_text(x = 100, y = 50, label = paste0("Km =", round(Km, 2), " ± ", round(Km.std, 2)))
  # If writing to file is desired 
  # ggsave(filename = paste0("/",filname,mm$enzymes[1],"_Michaelis_Menten.pdf"), plot = pl, device = "pdf")
  
  ll <- list(enzymes = mm$enzymes[1], 
             Km = Km, Km.std = Km.std, vmax = vmax, vmax.std = vmax.std,
             pl = pl, mml = mml)
  return(ll)
}

michaelisBoot <- function(mm, filname = NULL) {
  model.drm <- drm(v ~ S, data = mm, fct = MM.2())
  mml <- data.frame(S = seq(0, max(mm$S), length.out = 100))
  mml$v <- predict(model.drm, newdata = mml)
  
  # NLS

  bestfit <- nls(michaelis, mm, start=list(Vmax=6,Km=15))

  coef(bestfit)
  nlsb<-nlsBoot(bestfit, niter = 1000)

  plotfit(bestfit, smooth = TRUE)
  
  vmax <- nlsb$estiboot[1,1]
  vmax.std <- nlsb$estiboot[1,2]
  Km <- nlsb$estiboot[2,1]
  Km.std <- nlsb$estiboot[2,2]
  
  # Plot michaelis
  pl <- ggplot(mm, aes(x = S, y = v)) +
    theme_bw() +
    xlab("Substrate concentration") +
    ylab("Rate") +
    ggtitle(paste0("Michaelis-Menten kinetics of ", mm$enzymes[1])) +
    geom_point(alpha = 0.5) +
    geom_line(data = mml, aes(x = S, y = v), colour = "red") +
    geom_text(x = 50, y = 0.05, label = paste0("Vmax =", round(vmax, 2), " ± ", round(vmax.std, 2))) +
    geom_text(x = 50, y = 0.1, label = paste0("Km =", round(Km, 2), " ± ", round(Km.std, 2)))
  # If writing to file is desired
  #ggsave(filename = paste0("/", filname, mm$enzymes[1],"_Michaelis_Menten.pdf"), plot = pl, device = "pdf")
  
  # Plot from bootstrapping
  pdf(paste0("/",filname,mm$enzymes[1],"_plotfit.pdf"))
  plotfit(bestfit, smooth = TRUE)
  dev.off()
  
  # Plot michaelis
  pl <- ggplot(mm, aes(x = S, y = v)) +
    theme_bw() +
    xlab("Substrate concentration") +
    ylab("Rate") +
    ggtitle(paste0("Michaelis-Menten kinetics of ", mm$enzymes[1])) +
    geom_point() +
    geom_line(data = mml, aes(x = S, y = v), colour = "red") +
    geom_text(x = 50, y = 0.1, label = paste0("Vmax =", round(vmax, 2), " ± ", round(vmax.std, 2))) +
    geom_text(x = 50, y = 0.05, label = paste0("Km =", round(Km, 2), " ± ", round(Km.std, 2)))
  # If writing from file is desired
  # ggsave(filename = paste0("/",filname,mm$enzymes[1],"_Michaelis_Menten.pdf"), plot = pl, device = "pdf")
  
  ll <- list(enzymes = mm$enzymes[1], pl = pl, mml = mml, Km = Km, Km.std = Km.std, 
             vmax = vmax, vmax.std = vmax.std, boot = nlsb$estiboot)
  return(ll)
}
