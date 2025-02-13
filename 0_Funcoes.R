
# INTERPOLAÇÃO ------------------------------------------------------------

# Função para interpolar o hidrograma conforme um intervalo de tempo Δt especificado
fun.interpolar.hidrograma <- function(df.hidrograma,    # dataframe com duas colunas 't' e 'q'
                                      dt,               # intervalo de tempo Δt [min]
                                      metodo = "linear" # método de interpolação
                                      ){
  
  # Cria a função de interpolação
  interpolacao <- approxfun(x = df.hidrograma$t,
                            y = df.hidrograma$q,
                            method = metodo,
                            rule = 2) # mantém o último valor p/ extrapolação
  
  # Gera uma sequência de tempo discretizada em 'dt'
  t <- seq(from = 0, to = max(df.hidrograma$t), by = dt)
  
  # Garante que o último ponto (70 min) seja incluído
  if (max(t) < max(df.hidrograma$t)) {
    t <- c(t, max(df.hidrograma$t))
  }
  
  # Calcula vazões interpoladas
  q <- interpolacao(t)
  
  # Resultado
  return(data.frame(t = t, q = q))
  
}


# CHEZY -------------------------------------------------------------------

# Função objetivo de Chézy para ser minimizada pelo optimx
fun.chezy <- function(y, vazao, c.chezy, largura, decl){
  
  # Variáveis  
  q.hidro <- vazao  # vazão
  b <- largura      # largura do canal
  i.f <- decl       # declividade da linha de energia
  a <- b*y          # área molhada
  rh <- a/(2*y + b) # raio hidráulico
  
  # Vazão pela fórmula de chézy
  q.chezy <- c.chezy*a*sqrt(rh*i.f)
  
  # Definir função objetivo p/ encotrar a raiz
  fo <- q.hidro - q.chezy
  
  return(fo) 
  
}


# CARACTERÍSTICA NEGATIVA -------------------------------------------------

#  Função objetivo da característica negativa p/ ser minimzada e encontrar yp
fun.yp.montante <- function(par = c(NA),
                            vazao,
                            celeridade,
                            velocidade.s,
                            prof.s,
                            prof.m,
                            largura,
                            velocidade.m,
                            vazao.lateral,
                            decl.energia.m,
                            decl.canal,
                            dt){
  
  # Variáveis
  yp <- par[1]
  q <- vazao
  b <- largura
  cm <- celeridade
  vs <- velocidade.s
  vm <- velocidade.m
  ys <- prof.s
  ym <- prof.m
  ql <- vazao.lateral
  g <- 9.81
  ifm <- decl.energia.m
  i0 <- decl.canal
  am <- b*ym
  
  # Função de característica negativa c/ yp = q/(b*vp) substituída
  fo <- q/(yp*b) - vs - g/cm*(yp - ys) + ql/am*(vm + cm)*dt + g*(ifm - i0)*dt
  
  return(fo)
  
}

# Curva-chave
fun.yp.jusante <- function(par = c(NA),
                           celeridade,
                           velocidade.r,
                           prof.r,
                           prof.m,
                           largura,
                           velocidade.m,
                           vazao.lateral,
                           decl.energia.m,
                           decl.canal,
                           dt){
  
  # Variáveis
  yp <- par[1]
  b <- largura
  cm <- celeridade
  vr <- velocidade.r
  vm <- velocidade.m
  yr <- prof.r
  ym <- prof.m
  ql <- vazao.lateral
  g <- 9.81
  ifm <- decl.energia.m
  i0 <- decl.canal
  am <- b*ym
  
  # Função de característica positiva c/ curva chave subsituída em vp
  fo <- 111.7118*(yp - 0.3848)^1.2277/(b*yp) - vr + g/cm*(yp - yr) + ql/am*(vm - cm)*dt + g*(ifm - i0)*dt
  
  return(fo)
  
}
