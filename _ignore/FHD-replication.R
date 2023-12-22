
library(tidyverse)
library(MASS)
# library(Matrix)
library(expm)
options(scipen = 5) 
library(tidyverse)
library(knitr)
library(kableExtra)

ExR = 
    # Global Burden of Disease Collaborative Network. Global Burden of Disease Study 
    # 2019 (GBD 2019) Reference Life Table. Seattle, United States of America: 
    # Institute for Health Metrics and Evaluation (IHME), 2020
    tibble::tribble(
        ~Age, ~Life.Expectancy,
        0L,       88.8718951,
        1L,      88.00051053,
        5L,      84.03008056,
        10L,      79.04633476,
        15L,       74.0665492,
        20L,      69.10756792,
        25L,      64.14930031,
        30L,       59.1962771,
        35L,      54.25261364,
        40L,      49.31739311,
        45L,      44.43332057,
        50L,      39.63473787,
        55L,      34.91488095,
        60L,      30.25343822,
        65L,      25.68089534,
        70L,      21.28820012,
        75L,      17.10351469,
        80L,      13.23872477,
        85L,      9.990181244,
        90L,      7.617724915,
        95L,      5.922359078
    )

fExR <- function(x) pmax(0,unname(Hmisc::approxExtrap(ExR$Age, ExR$Life.Expectancy,xout = x)$y))

alt_simp_coef <- function(i) c(17, 59, 43, 49, rep(48, i-8), 49, 43, 59, 17) / 48
alt_simp      <- function(x,h) h*sum(alt_simp_coef(length(x)) * x)

age0 =35
horizon = 10
r_ann = .03
Delta_t = 1/12
r_Delta_t = r_ann * Delta_t #((1 + r_ann)^Delta_t) - 1
cycles = horizon/Delta_t

#FRH consider a woman who develops bipolar depression at age 35, 
#lives for 10 years with the disorder, and then dies prematurely 
#at age 45. The disability weight is 0.60, life expectancy at age 45 
#is 34.73 years (age 79.73), and the “present” for the present value 
#calculation is age 35 (meaning exactly on her 35th birthday).

Q = list()
for (i in 1:(cycles)) {
    Q = append(Q,{
        list(matrix(c(0,0,0,0,
                 0,0,0,0,
                 0,0,0,0,
                 0,0,0,0
        ),
        nrow = 4, 
        ncol=4,
        byrow=TRUE,
        dimnames = list(c("H","S","D","trD"),c("H","S","D","trD"))))
    })
}

Q = append(Q,list(matrix(c(0,0,0,0,
                  0,-100,100,100,
                  0,0,0,0,
                  0,0,0,0
),
nrow = 4, 
ncol=4,
byrow=TRUE,
dimnames = list(c("H","S","D","trD"),c("H","S","D","trD")))))

P = Q %>% map(~({expm(.x)}))
P = P %>% map(~({
    tmp_ = .x
    tmp_["trD","trD"] = 0
    tmp_
}))

s0T = c(0,1,0,0) 

s = s0T 
res = matrix(0,nrow=11,ncol=4)
trace = 
    P %>% map_df(~({
    s <<- s %*% .x
    data.frame(s)
})) 

# Payoff Vectors
d = c(0,0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,0); d
#d = c(0,.6,0,0)
et = 0:cycles %>% map(~({
    if (.x==cycles) c(0,0,0,(1/r_ann)*(1-exp(-r_ann*34.73))) else c(0,0,0,0)
    }))
y = 0:cycles %>% map(~({
    if (.x==cycles) c(0,Delta_t * 0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,(1/r_ann)*(1-exp(-r_ann*34.73))) else c(0,Delta_t * 0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,0)
}))

YLDt = as.vector(as.matrix(trace) %*% d) * Delta_t 
YLLt = map2_dbl(0:cycles,et,~({as.matrix(trace[.x+1,]) %*% as.matrix(.y)}))
DALYt = map2_dbl(0:cycles,y,~({as.matrix(trace[.x+1,]) %*% as.matrix(.y)}))
disc_discrete = 0:cycles %>% map_dbl(~({1/(1+r_Delta_t)^.x}))
disc_cont = 0:cycles %>% map_dbl(~({exp(-r_ann * Delta_t * .x)}))

trace_final <- 
    trace %>% 
    as.data.frame() %>% 
    mutate(YLDt = YLDt,
           dYLDt = disc_cont * YLDt,
           ddYLDt = disc_discrete * YLDt,
           YLLt = YLLt, 
           dYLLt = disc_cont * YLLt,
           ddYLLt = disc_discrete * YLLt,
           DALYt = DALYt,
           dDALYt  = disc_cont * DALYt,
           ddDALYt = disc_discrete * DALYt) %>% 
    mutate(DALYt_ = YLDt + YLLt,
           dDALYt_ = dYLDt + dYLLt,
           ddDALYt_ = ddYLDt + ddYLLt) %>% 
    mutate(disc_cont = disc_cont,
           disc_discrete = disc_discrete) %>% 
    mutate(at = (row_number()-1)*Delta_t+age0) 
    
YLD_FRH = (1/r_ann)*0.6*(1-exp(-r_ann*10))
YLL_FRH = exp(-r_ann*10)*(1/r_ann)*(1-exp(-r_ann*34.73))

trace_final %>% 
    summarise_at(vars(YLDt,dYLDt,ddYLDt),sum)  %>% 
    mutate(FRH = YLD_FRH)

trace_final %>% 
    summarise_at(vars(YLLt,dYLLt,ddYLLt),sum)  %>% 
    mutate(FRH = YLL_FRH)

trace_final %>% 
    summarise_at(vars(DALYt,dDALYt,ddDALYt),sum)  %>% 
    mutate(FRH = YLD_FRH + YLL_FRH)
    



### APPROACH 2

# Payoff Vectors
d = c(0,0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,0); d
et = 0:cycles %>% map(~({
    if (.x==cycles) c(0,0,0,(1/r_ann)*(1-exp(-r_ann*34.73))) else c(0,0,0,0)
}))
y = 0:cycles %>% map(~({
    if (.x==cycles) c(0,Delta_t * 0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,(1/r_ann)*(1-exp(-r_ann*34.73))) else c(0,Delta_t * 0.6 * (1/r_Delta_t)*(1-exp(-r_Delta_t)),0,0)
}))

YLDt = as.vector(as.matrix(trace) %*% d) * Delta_t 
YLLt = map2_dbl(0:cycles,et,~({as.matrix(trace[.x+1,]) %*% as.matrix(.y)}))
DALYt = map2_dbl(0:cycles,y,~({as.matrix(trace[.x+1,]) %*% as.matrix(.y)}))
disc_discrete = 0:cycles %>% map_dbl(~({1/(1+r_Delta_t)^.x}))
disc_cont = 0:cycles %>% map_dbl(~({exp(-r_ann * Delta_t * .x)}))

trace_final <- 
    trace %>% 
    as.data.frame() %>% 
    mutate(YLDt = YLDt,
           dYLDt = disc_cont * YLDt,
           ddYLDt = disc_discrete * YLDt,
           YLLt = YLLt, 
           dYLLt = disc_cont * YLLt,
           ddYLLt = disc_discrete * YLLt,
           DALYt = DALYt,
           dDALYt  = disc_cont * DALYt,
           ddDALYt = disc_discrete * DALYt) %>% 
    mutate(DALYt_ = YLDt + YLLt,
           dDALYt_ = dYLDt + dYLLt,
           ddDALYt_ = ddYLDt + ddYLLt) %>% 
    mutate(disc_cont = disc_cont,
           disc_discrete = disc_discrete) %>% 
    mutate(at = (row_number()-1)*Delta_t+age0) 

YLD_FRH = (1/r_ann)*0.6*(1-exp(-r_ann*10))
YLL_FRH = exp(-r_ann*10)*(1/r_ann)*(1-exp(-r_ann*34.73))

trace_final %>% 
    summarise_at(vars(YLDt,dYLDt,ddYLDt),sum)  %>% 
    mutate(FRH = YLD_FRH)


trace_final %>% 
    summarise_at(vars(YLLt,dYLLt,ddYLLt),sum)  %>% 
    mutate(FRH = YLL_FRH)

trace_final %>% 
    summarise_at(vars(DALYt,dDALYt,ddDALYt),sum)  %>% 
    mutate(FRH = YLD_FRH + YLL_FRH)


    


