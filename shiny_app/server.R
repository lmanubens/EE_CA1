### EE CA1 Shiny app 
  
  library(shiny)
  library(plotly)
  library(ggpubr)
  library(afex)
  library(dplyr)
  library(reshape2)
  library(xtable)
  library(emmeans)
  # library(Rmpfr)
  library(ggsignif)
  library(plotrix)
  library(scales)
  # library(forcats)
  # library(Hmisc)
  library(rstatix)
  
  shinyServer(function(input, output, session) {
    # upDataorig <- reactive({
    #   load('subsetdata.Rdata')
    #   return(data)
    # })
    # upDatagroups <- reactive({
    #   groupsdf <- read.csv('freqinout_rad.csv', header = TRUE)
    #   return(groupsdf)
    # })
    observe({
      # updateSelectInput(
      #   session,
      #   "variablegroups",
      #   choices=names(upDatagroups()),
      #   selected=names(upDatagroups())[1]
      # )
      ##############################
      # if(input$selectall == 0) {
        updateCheckboxGroupInput(
          session,
          "variablegroups",
          choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
          selected=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE')
          )
        updateSelectInput(
          session,
          "variableshollgroup",
          choices=c('EE 2mo', 'EE 6mo', 'Age effect'),
          selected='EE 2mo'
        )
        updateSelectInput(
          session,
          "variableinoutgroup",
          choices=c('EE 2mo', 'EE 6mo', 'Age effect'),
          selected='EE 2mo'
        )
        
        updateSelectInput(
          session,
          "variableinoutregion",
          choices=c('Str. Radiatum', 'Both', 'Str. Radiatum - only morph.','Both - only morph.'),
          selected='Str. Radiatum'
        )
        updateSelectInput(
          session,
          "variablecomp",
          choices=c('EE 2mo','EE 6mo','Age effect'),
          selected='EE 2mo'
        )
        updateSelectInput(
          session,
          "variablecomp2",
          choices=c('EE 2mo','EE 6mo','Age effect'),
          selected='EE 2mo'
        )
        updateSelectInput(
          session,
          "variablecomp3",
          choices=c('EE 2mo','EE 6mo','Age effect'),
          selected='EE 2mo'
        )
        updateSelectInput(
          session,
          "variablecomp4",
          choices=c('EE 2mo','EE 6mo','Age effect'),
          selected='EE 2mo'
        )
        updateSelectInput(
          session,
          "variablemet",
          choices=c('DYRK1A_act','OR','spine_tot','spine_A','spine_B','spine_C','spine_B&C','vglutvgat','ratio_fEPSP'),
          selected='DYRK1A_act'
        )
        updateSelectInput(
          session,
          "variablerecmet",
          choices=c('tot_len','max_plen','bpoints','mpeucl','maxbo','mangleB','mblen',
                    'mplen','mbo','wh','wz','chullx','chully','chullz','hull','masym',
                    'mparea','span_vol','DV_span','AUC','con_rep','con_rep_onlymorph',
                    'Succ_spikes_rad_40Hz','Succ_spikes_both_40Hz','Succ_spikes_rad_100Hz','Succ_spikes_both_100Hz'),
          selected='DV_span'
        )
        updateSelectInput(
          session,
          "variablesubl",
          choices=c('Apical CA1','Stratum Radiatum','Stratum Lacunosum'),
          selected='Apical CA1'
        )
        # updateSelectInput(
        #   session,
        #   "variablegroupsdisc",
        #   choices=c('Acute EE (2mo)','Discontinued EE (6mo)'),
        #   selected= 'Acute EE (2mo)'
        # )
        
      # }
      # else if (input$selectall%%2 == 0)
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variablegroups",
      #     choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
      #     selected=c()
      #   )
      # }
      # else
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variablegroups",
      #     choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
      #     selected=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE')
      #   )
      # }
      
      # if(input$selectall2 == 0) {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variablegroups2",
      #     choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
      #     selected=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE')
      #   )
      # }
      # else if (input$selectall2%%2 == 0)
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variablegroups2",
      #     choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
      #     selected=c()
      #   )
      # }
      # else
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variablegroups2",
      #     choices=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE'),
      #     selected=c('WT','TgDyrk1A','WT-EE','TgDyrk1A-EE','WT-NE_NE','TgDyrk1A-NE_NE','WT-EE_NE','TgDyrk1A-EE_NE')
      #   )
      # }
      # ############################################################
      
      ##############################################################
      # if(input$selectall2 == 0) {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variabledist",
      #     choices=names(upDataorig()[22:28]),
      #     selected=names(upDataorig())[c(22,28)]
      #   )
      # }
      # else if (input$selectall2%%2 == 0)
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variabledist",
      #     choices=names(upDataorig()[22:28]),
      #     selected=c()
      #   )
      # }
      # else 
      # {
      #   updateCheckboxGroupInput(
      #     session,
      #     "variabledist",
      #     choices=names(upDataorig()[22:28]),
      #     selected=names(upDataorig()[22:28])
      #   )
      # }
      
    })
    # upData <- reactive({
    #   
    #   dat <- subset(upDataorig(),select = c(input$variablemorph,input$variabledist,input$variableiq))
    #   return(dat)
    # })
    ###################
    ###################
    # output$freqinout <- renderPlotly({
    #   if(input$variableinoutgroup=='Proximal - Only morphology'){
    #     df <- read.csv('freqinout_onlymorph_rad_prox.csv', header = TRUE)
    #   }
    #   else if(input$variableinoutgroup=='Distal - Only morphology'){
    #     df <- read.csv('freqinout_onlymorph_rad_dist.csv', header = TRUE)
    #   }
    #   else if(input$variableinoutgroup=='All Stratum Radiatum - Only morphology'){
    #     df <- read.csv('freqinout_onlymorph_rad.csv', header = TRUE)
    #   }
    #   else if(input$variableinoutgroup=='Proximal'){
    #     df <- read.csv('freqinout_rad_prox.csv', header = TRUE)
    #   }
    #   else if(input$variableinoutgroup=='Distal'){
    #     df <- read.csv('freqinout_rad_dist.csv', header = TRUE)
    #   }
    #   else if(input$variableinoutgroup=='All Stratum Radiatum'){
    #     df <- read.csv('freqinout_rad.csv', header = TRUE)
    #   }
    #   df$genotype <- rep(c(rep('Wild-Type',18), #B6EiC3Sn
    #                        rep('TgDyrk1A',19),
    #                        rep('Wild-Type_EE',19),
    #                        rep('TgDyrk1A_EE',15),
    #                        rep('Wild-Type_NE_NE',18),
    #                        rep('TgDyrk1A_NE_NE',18),
    #                        rep('Wild-Type_EE_NE',18),
    #                        rep('TgDyrk1A_EE_NE',18)),10)
    #   df$genotype2 <- rep(c(rep('WT',18), #B6EiC3Sn
    #                        rep('TgDyrk1A',19),
    #                        rep('WT-EE',19),
    #                        rep('TgDyrk1A-EE',15),
    #                        rep('WT-NE_NE',18),
    #                        rep('TgDyrk1A-NE_NE',18),
    #                        rep('WT-EE_NE',18),
    #                        rep('TgDyrk1A-EE_NE',18)),10)
    #   
    #   df$genotype <- factor(df$genotype,
    #                             levels=c('Wild-Type',
    #                                      'TgDyrk1A',
    #                                      'Wild-Type_EE',
    #                                      'TgDyrk1A_EE',
    #                                      'Wild-Type_NE_NE',
    #                                      'TgDyrk1A_NE_NE',
    #                                      'Wild-Type_EE_NE',
    #                                      'TgDyrk1A_EE_NE'))
    #   
    #   dfplot <- filter(df, genotype2 %in% input$variablegroups)
    #   
    #   
    #   p <- ggline(dfplot, x = "FreqIn", y = "FreqOut", add = c("mean_se"),#,"point"),
    #               color = "genotype", 
    #               size=1.5,
    #               add.params = list(alpha=0.5, size=0.7),
    #               shape = "genotype",
    #               linetype = "genotype",
    #               palette = "Set1",
    #               ggtheme = theme_bw(),
    #               # title='Input Output Relation',
    #               xlab='Input frequency [Hz]',
    #               ylab='Output frequency [Hz]',
    #               xlim=c(0,100),
    #               ylim=c(0,100),
    #               legend='right') +
    #     theme(axis.text=element_text(size=20),
    #           axis.title = element_text(size=25),
    #           title=element_text(size=25),
    #           legend.text = element_text(size=25)) +
    #     stat_compare_means(aes(group = genotype),
    #                        method="kruskal.test",
    #                        label = "p.signif",
    #                        label.y=100,
    #                        hide.ns = T,
    #                        size=8)
    #   p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
    #   p
    #   # print(
    #   #   ggplotly(p)
    #   # )
    # })
    
    ###################
    ###################
    
    plotfreqinout <- reactive({
      
      if(input$variableinoutregion == "Str. Radiatum - only morph."){
        df <- read.csv('freqinout_rad_BC_onlymoprh.csv', header = TRUE)
      }
      else if(input$variableinoutregion == "Str. Radiatum"){
        df <- read.csv('freqinout_rad_BC.csv', header = TRUE)
      }
      else if(input$variableinoutregion == "Both - only morph."){
        df <- read.csv('freqinout_both_BC_onlymorph.csv', header = TRUE)
      }
      else if(input$variableinoutregion == "Both"){
        df <- read.csv('freqinout_both_BC.csv', header = TRUE)
      }
      
      df <- df[,1:3]
      df$genotype <- as.factor(rep(c(rep('WTNE',18), #B6EiC3Sn
                           rep('TGNE',19),
                           rep('WTEE',19),
                           rep('TGEE',15),
                           rep('WTNE_NE',18),
                           rep('TGNE_NE',18),
                           rep('WTEE_NE',18),
                           rep('TGEE_NE',18)),10))
      blocking <- read.csv('blocking.csv',header=F)
      df$blocking <- as.factor(rep(blocking$V1,10))
      
      if(input$variableinoutgroup == "EE 2mo"){
        df <- subset(df,df$genotype=='WTNE'|df$genotype=='WTEE'|df$genotype=='TGNE'|df$genotype=='TGEE')
      }
      else if(input$variableinoutgroup == "EE 6mo"){
        df <- subset(df,df$genotype=='WTNE_NE'|df$genotype=='WTEE_NE'|df$genotype=='TGNE_NE'|df$genotype=='TGEE_NE')
        levels(df$genotype)[levels(df$genotype)=="WTNE_NE"] <- "WTNE"
        levels(df$genotype)[levels(df$genotype)=="TGNE_NE"] <- "TGNE"
        levels(df$genotype)[levels(df$genotype)=="WTEE_NE"] <- "WTEE"
        levels(df$genotype)[levels(df$genotype)=="TGEE_NE"] <- "TGEE"
      }
      else{
        dftot <- df
        df <- subset(df,df$genotype=='WTNE'|df$genotype=='TGNE'|df$genotype=='WTEE'|df$genotype=='TGEE')
        df$gen <- factor(c('WT','TG'))
        df$age <- NULL
        df$treat <- NULL
        df$gen[df$genotype=="WTNE"] <- 'WT'
        df$gen[df$genotype=="WTEE"] <- 'WT'
        df$gen[df$genotype=="TGNE"] <- 'TG'
        df$gen[df$genotype=="TGEE"] <- 'TG'
        df$treat[df$genotype=="WTNE"] <- 'NE'
        df$treat[df$genotype=="WTEE"] <- 'EE'
        df$treat[df$genotype=="TGNE"] <- 'NE'
        df$treat[df$genotype=="TGEE"] <- 'EE'
        df$age <- '2mo'
        df$genotype <- as.factor(df$genotype)
        # df <- subset(df,df$genotype=='WTNE' | df$genotype=='TGNE')
        levels(df$genotype)[levels(df$genotype)=="WTNE"] <- 'WTNE_2mo'
        levels(df$genotype)[levels(df$genotype)=="TGNE"] <- 'TGNE_2mo'
        # levels(df$genotype)[levels(df$genotype)=="WTEE"] <- 'WTEE'
        # levels(df$genotype)[levels(df$genotype)=="TGEE"] <- 'TGEE'
        df1 <- df
        
        df <- subset(dftot,dftot$genotype=='WTNE_NE'|dftot$genotype=='TGNE_NE'|dftot$genotype=='WTEE_NE'|dftot$genotype=='TGEE_NE')
        df$gen <- factor(c('WT','TG'))
        df$age <- NULL
        df$treat <- NULL
        df$gen[df$genotype=="WTNE_NE"] <- 'WT'
        df$gen[df$genotype=="WTEE_NE"] <- 'WT'
        df$gen[df$genotype=="TGNE_NE"] <- 'TG'
        df$gen[df$genotype=="TGEE_NE"] <- 'TG'
        df$treat[df$genotype=="WTNE_NE"] <- 'NE'
        df$treat[df$genotype=="WTEE_NE"] <- 'EE'
        df$treat[df$genotype=="TGNE_NE"] <- 'NE'
        df$treat[df$genotype=="TGEE_NE"] <- 'EE'
        df$age <- '6mo'
        df$genotype <- as.factor(df$genotype)
        # df <- subset(df,df$genotype=='WTNE' | df$genotype=='TGNE') 
        levels(df$genotype)[levels(df$genotype)=="WTNE_NE"] <- 'WTNE_6mo'
        levels(df$genotype)[levels(df$genotype)=="TGNE_NE"] <- 'TGNE_6mo'
        levels(df$genotype)[levels(df$genotype)=="WTEE_NE"] <- 'WTEEdis'
        levels(df$genotype)[levels(df$genotype)=="TGEE_NE"] <- 'TGEEdis'
        
        # levels(df$genotype)[levels(df$genotype)=="WTNE"] <- "WT2mo"
        # levels(df$genotype)[levels(df$genotype)=="TGNE"] <- "TG2mo"
        # levels(df$genotype)[levels(df$genotype)=="WTNE_NE"] <- "WT6mo"
        # levels(df$genotype)[levels(df$genotype)=="TGNE_NE"] <- "TG6mo"
        df <- rbind(df1,df)
        
        df$genotype <- factor(df$genotype,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
      }
      if(input$variableinoutgroup != "Age effect"){
        df$genotype <- factor(df$genotype,
                              levels=c('WTNE',
                                       'TGNE',
                                       'WTEE',
                                       'TGEE'))
      }else{
        # df$genotype <- factor(df$genotype,
        #                       levels=c('WT2mo',
        #                                'WT6mo',
        #                                'TG2mo',
        #                                'TG6mo'))
      }
      
      if(input$variableinoutgroup != "Age effect"){
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$treat <- sapply(df$genotype,substring,first=3,last=4)
        
        annogen <- NULL
        annotreat <- NULL
        annogentreat <- NULL
        annowttg <- NULL
        annowtwtee <- NULL
        annowttgee <- NULL
        annotgtgee <- NULL
        # for(i in seq(10,100,10)){
        #   dfi <- subset(df,df$FreqIn==i)
        #   print(summary(mixed(FreqOut~gen*treat+(1|blocking),data=dfi)))
        #   groupannop <- as.data.frame(summary(mixed(FreqOut~gen*treat+(1|blocking),data=dfi))[[10]])[2:4,5]
        #   annogen[i/10] <- groupannop[1]
        #   annotreat[i/10] <- groupannop[2]
        #   annogentreat[i/10] <- groupannop[3]
        #   annowttg[i/10] <- as.data.frame(pairs(emmeans(mixed(FreqOut~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[6]
        #   annowtwtee[i/10] <- as.data.frame(pairs(emmeans(mixed(FreqOut~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[5]
        #   annowttgee[i/10] <- as.data.frame(pairs(emmeans(mixed(FreqOut~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[3]
        #   annotgtgee[i/10] <- as.data.frame(pairs(emmeans(mixed(FreqOut~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[2]
        # }
        # 
        # stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[1:10]), corr = FALSE,
        #                                                              cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                              symbols = c("****", "***", "**", "*", "")))),
        #                             ptreat=c('Treat',as.character(symnum(as.numeric(annotreat[1:10]), corr = FALSE,
        #                                                                  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                  symbols = c("****", "***", "**", "*", "")))),
        #                             pgentreat=c('Gen:Treat',as.character(symnum(as.numeric(annogentreat[1:10]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pwttg=c('WT-TG',as.character(symnum(as.numeric(annowttg[1:10]), corr = FALSE,
        #                                                                 cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                 symbols = c("****", "***", "**", "*", "")))),
        #                             pwtwtee=c('WT-WTEE',as.character(symnum(as.numeric(annowtwtee[1:10]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))),
        #                             pwttgee=c('WT-TGEE',as.character(symnum(as.numeric(annowttgee[1:10]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))),
        #                             ptgtgee=c('TG-TGEE',as.character(symnum(as.numeric(annotgtgee[1:10]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))))
      }else{
        # df$gen <- sapply(df$genotype,substring,first=1,last=2)
        # df$age <- sapply(df$genotype,substring,first=3,last=5)
        
        annogen <- NULL
        annoage <- NULL
        annogenage <- NULL
        annoagetreat <- NULL
        annogenagetreat <- NULL
        annodisc <- NULL
        annowtage <- NULL
        annotgage <- NULL
        annowteedis <- NULL
        annotgeedis <- NULL
        
        # for(i in seq(10,100,10)){
        #   dfi <- subset(df,df$FreqIn==i)
        #   tryCatch(print(summary(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi))), error=function(e){})
        #   groupannop <- tryCatch(as.data.frame(summary(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi))[[10]])[2:8,5],error=function(e){c(1,1,1)})
        #   annogen[i/10] <- groupannop[1]
        #   agenep <- tryCatch(as.data.frame(joint_tests(mixed(value~gen*age*treat+(1|blocking),data=dfi), by = "treat"))$p.value[4],error=function(e){1})
        #   annoage[i/10] <- agenep  #groupannop[2]
        #   annogenage[i/10] <- groupannop[3]
        #   annoagetreat[i/10] <- groupannop[6]
        #   annogenagetreat[i/10] <- groupannop[7]
        #   discannop <- tryCatch(as.data.frame(joint_tests(mixed(value~gen*age*treat+(1|blocking),data=dfi), by = "treat"))$p.value[3],error=function(e){1})
        #   annodisc[i/10] <- discannop 
        #   annowtage[i/10] <- tryCatch(as.data.frame(pairs(emmeans(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[27],error=function(e){1})
        #   annotgage[i/10] <- tryCatch(as.data.frame(pairs(emmeans(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[24],error=function(e){1})
        #   annowteedis[i/10] <- tryCatch(as.data.frame(pairs(emmeans(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[9],error=function(e){1})
        #   annotgeedis[i/10] <- tryCatch(as.data.frame(pairs(emmeans(mixed(FreqOut~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[2],error=function(e){1})
        # }
        
        # stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[1:10]), corr = FALSE,
        #                                                              cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                              symbols = c("****", "***", "**", "*", "")))),
        #                             page=c('Age NE',as.character(symnum(as.numeric(annoage[1:10]), corr = FALSE,
        #                                                              cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                              symbols = c("****", "***", "**", "*", "")))),
        #                             pgenage=c('Gen:Age',as.character(symnum(as.numeric(annogenage[1:10]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))),
        #                             pagetreat=c('Disc',as.character(symnum(as.numeric(annodisc[1:10]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pgenagetreat=c('Gen:Age:Treat',as.character(symnum(as.numeric(annogenagetreat[1:10]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pwtage=c('WT-age',as.character(symnum(as.numeric(annowtage[1:10]), corr = FALSE,
        #                                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                   symbols = c("****", "***", "**", "*", "")))),
        #                             ptgage=c('TG-age',as.character(symnum(as.numeric(annotgage[1:10]), corr = FALSE,
        #                                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                   symbols = c("****", "***", "**", "*", "")))),
        #                             pwteedis=c('WT-EEdis',as.character(symnum(as.numeric(annowteedis[1:10]), corr = FALSE,
        #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                       symbols = c("****", "***", "**", "*", "")))),
        #                             ptgeedis=c('TG-EEdis',as.character(symnum(as.numeric(annotgeedis[1:10]), corr = FALSE,
        #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                       symbols = c("****", "***", "**", "*", "")))))
      }
      # df <- subset(df,df$time!=0)
      # stat.pairwise <- stat.pairwise[1:11,]
      
      
      if(input$variableinoutgroup != "Age effect"){
        p <- ggline(df, x = "FreqIn", y = "FreqOut", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=1.5,
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    palette = "Set1",
                    ggtheme = theme_bw(),
                    # title='Input Output Relation',
                    xlab='Input frequency [Hz]',
                    ylab='Output frequency [Hz]',
                    xlim=c(0,100),
                    ylim=c(0,120),
                    legend='right') +
          theme(axis.text=element_text(size=20),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25)) 
        p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
        
        p #+ 
          # geom_signif(annotations=stat.pairwise$pgentreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 115,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pgen,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 110,size=0,color='orange3') +
          # geom_signif(annotations=stat.pairwise$ptreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 105,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwttg,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 95,size=0,color='#377eb8') +
          # geom_signif(annotations=stat.pairwise$pwtwtee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 90,size=0,color='#4daf4a') +
          # geom_signif(annotations=stat.pairwise$pwttgee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 80,size=0,color='#f781bf') +
          # geom_signif(annotations=stat.pairwise$ptgtgee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 75,size=0,color='#984ea3')  
      }else{
        p <- ggline(df, x = "FreqIn", y = "FreqOut", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=1.5,
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    palette = "Set2",
                    ggtheme = theme_bw(),
                    # title='Input Output Relation',
                    xlab='Input frequency [Hz]',
                    ylab='Output frequency [Hz]',
                    xlim=c(0,100),
                    ylim=c(0,105),
                    legend='right') +
          theme(axis.text=element_text(size=20),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25)) 
        p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
        p #+ 
          # geom_signif(annotations=stat.pairwise$pagetreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 100,size=0,color='brown') +
          # # geom_signif(annotations=stat.pairwise$pgenage,textsize=4,aes(group = genotype),
          # #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          # #             y_position = 95,size=0,color='brown') +
          # # geom_signif(annotations=stat.pairwise$pgenagetreat,textsize=4,aes(group = genotype),
          # #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          # #             y_position = 90,size=0,color='orange3') +
          # geom_signif(annotations=stat.pairwise$page,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 90,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwteedis,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 80,size=0,color='#66c2a5') +
          # geom_signif(annotations=stat.pairwise$ptgeedis,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 85,size=0,color='#8da0cb') +
          # geom_signif(annotations=stat.pairwise$pwtage,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 70,size=0,color='#66c2a5') +
          # geom_signif(annotations=stat.pairwise$ptgage,textsize=4,aes(group = genotype),
          #             xmin=seq(0,100,10),xmax=seq(0,100,10),tip_length=0,
          #             y_position = 75,size=0,color='#8da0cb')
      }
      
      # print(
      #   ggplotly(p)
      # )
    })
    
    output$freqinout <-  renderPlot({
      print(plotfreqinout())
    })
    
    
    output$downloadPlot5 <- downloadHandler(
      filename = function(){paste('FreqInOut_',input$variableinoutgroup, '.pdf', sep = '')},
      
      content = function(file){
        ggsave(file,width=14,height=12,plotfreqinout(),useDingbats=F)
      },
      
      contentType = "application/pdf"
    )
    
    ###################
    ###################
    plotslopeltp <- reactive({
      if(input$variablecomp4 == "EE 2mo"){
        df <- read.csv('slopes_2mo.csv', header = TRUE)
      }
      else if(input$variablecomp4 == "EE 6mo"){
        df <- read.csv('slopes_6mo.csv', header = TRUE)
      }
      else{
        data1 <- read.csv('slopes_2mo.csv', header = TRUE)
        # df1 <- subset(df1,df1$genotype=='WTNE' | df1$genotype=='TGNE')
        # levels(df1$genotype)[levels(df1$genotype)=="WTNE"] <- "WT2mo"
        # levels(df1$genotype)[levels(df1$genotype)=="TGNE"] <- "TG2mo"
        data2 <- read.csv('slopes_6mo.csv', header = TRUE)
        # df2 <- subset(df2,df2$genotype=='WTNE' | df2$genotype=='TGNE')
        # levels(df2$genotype)[levels(df2$genotype)=="WTNE"] <- "WT6mo"
        # levels(df2$genotype)[levels(df2$genotype)=="TGNE"] <- "TG6mo"
        # df <- rbind(df1,df2)
        
        data1$gen <- NULL
        data1$age <- NULL
        data1$treat <- NULL
        data1$gen[data1$genotype=="WTNE"] <- 'WT'
        data1$gen[data1$genotype=="WTEE"] <- 'WT'
        data1$gen[data1$genotype=="TGNE"] <- 'TG'
        data1$gen[data1$genotype=="TGEE"] <- 'TG'
        data1$treat[data1$genotype=="WTNE"] <- 'NE'
        data1$treat[data1$genotype=="WTEE"] <- 'EE'
        data1$treat[data1$genotype=="TGNE"] <- 'NE'
        data1$treat[data1$genotype=="TGEE"] <- 'EE'
        data1$age <- '2mo'
        
        data1$genotype <- as.factor(data1$genotype)
        # data1 <- subset(data1,data1$genotype=='WTNE' | data1$genotype=='TGNE') 
        levels(data1$genotype)[levels(data1$genotype)=="WTNE"] <- 'WTNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="TGNE"] <- 'TGNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="WTEE"] <- 'WTEE'
        levels(data1$genotype)[levels(data1$genotype)=="TGEE"] <- 'TGEE'
        
        data2$gen <- NULL
        data2$age <- NULL
        data2$treat <- NULL
        data2$gen[data2$genotype=="WTNE"] <- 'WT'
        data2$gen[data2$genotype=="WTEE"] <- 'WT'
        data2$gen[data2$genotype=="TGNE"] <- 'TG'
        data2$gen[data2$genotype=="TGEE"] <- 'TG'
        data2$treat[data2$genotype=="WTNE"] <- 'NE'
        data2$treat[data2$genotype=="WTEE"] <- 'EE'
        data2$treat[data2$genotype=="TGNE"] <- 'NE'
        data2$treat[data2$genotype=="TGEE"] <- 'EE'
        data2$age <- '6mo'
        data2$genotype <- as.factor(data2$genotype)
        # data2 <- subset(data2,data2$genotype=='WTNE' | data2$genotype=='TGNE') 
        levels(data2$genotype)[levels(data2$genotype)=="WTNE"] <- 'WTNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="TGNE"] <- 'TGNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="WTEE"] <- 'WTEEdis'
        levels(data2$genotype)[levels(data2$genotype)=="TGEE"] <- 'TGEEdis'
        
        df <- rbind(data1,data2)
        
        # data$genotype <- factor(data$genotype,levels=c('Wild-Type_2mo','Wild-Type_6mo','TgDyrk1A_2mo','TgDyrk1A_6mo'))
        df$genotype <- factor(df$genotype,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
      }
      if(input$variablecomp4 != "Age effect"){
        df$genotype <- factor(df$genotype,
                              levels=c('WTNE',
                                       'TGNE',
                                       'WTEE',
                                       'TGEE'))
      }else{
        # df$genotype <- factor(df$genotype,
        #                       levels=c('WT2mo',
        #                                'WT6mo',
        #                                'TG2mo',
        #                                'TG6mo'))
      }
      
      df <- filter(df, df$time %in% seq(-15,60,1) )
      # df[is.na(df$slope),]$slope <- 0
      df <- remove_missing(df)
      
      if(input$variablecomp4 != "Age effect"){
        # df$gen <- sapply(df$genotype,substring,first=1,last=2)
        # df$treat <- sapply(df$genotype,substring,first=3,last=4)
        
        annogen <- NULL
        annotreat <- NULL
        annogentreat <- NULL
        annowttg <- NULL
        annowtwtee <- NULL
        annowttgee <- NULL
        annotgtgee <- NULL
        # for(i in seq(-15,60,1)){
        #   dfi <- subset(df,df$time==i)
        #   print(summary(aov(slope~gen*treat,data=dfi)))
        #   groupannop <- as.data.frame(summary(aov(slope~gen*treat,data=dfi))[[1]])[1:3,5]
        #   annogen[i+16] <- groupannop[1]
        #   annotreat[i+16] <- groupannop[2]
        #   annogentreat[i+16] <- groupannop[3]
        #   annowttg[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[6]
        #   annowtwtee[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[5]
        #   annowttgee[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[3]
        #   annotgtgee[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[2]
        # }
        
        # stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[2:76]), corr = FALSE,
        #                                               cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                               symbols = c("****", "***", "**", "*", "")))),
        #                             ptreat=c('Treat',as.character(symnum(as.numeric(annotreat[2:76]), corr = FALSE,
        #                                                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                      symbols = c("****", "***", "**", "*", "")))),
        #                             pgentreat=c('Gen:Treat',as.character(symnum(as.numeric(annogentreat[2:76]), corr = FALSE,
        #                                                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                      symbols = c("****", "***", "**", "*", "")))),
        #                             pwttg=c('WT-TG',as.character(symnum(as.numeric(annowttg[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pwtwtee=c('WT-WTEE',as.character(symnum(as.numeric(annowtwtee[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pwttgee=c('WT-TGEE',as.character(symnum(as.numeric(annowttgee[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             ptgtgee=c('TG-TGEE',as.character(symnum(as.numeric(annotgtgee[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))))
      }else{
        # df$gen <- sapply(df$genotype,substring,first=1,last=2)
        # df$age <- sapply(df$genotype,substring,first=3,last=5)
        
        annogen <- NULL
        annoage <- NULL
        annogenage <- NULL
        annoagetreat <- NULL
        annogenagetreat <- NULL
        annowtage <- NULL
        annotgage <- NULL
        annowteedis <- NULL
        annotgeedis <- NULL
       #  for(i in seq(-15,60,1)){
       #    dfi <- subset(df,df$time==i)
       #    # print(summary(anova_test(slope~gen*age*treat,data=dfi)))
       #    print(anova_test(slope~gen*age*treat,data=dfi))
       #    # groupannop <- as.data.frame(anova_test(slope~gen*age*treat,data=dfi)[[1]])[1:3,5]
       #    groupannop <- as.data.frame(anova_test(slope~gen*age*treat,data=dfi))$p
       #    ################### HERE tbc
       #    annogen[i+16] <- groupannop[1]
       #    annoage[i+16] <- groupannop[2]
       #    annogenage[i+16] <- groupannop[4]
       #    annoagetreat[i+16] <- groupannop[6]
       #    annogenagetreat[i+16] <- groupannop[7]
       #    annowtage[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*age*treat,data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[27]
       #    annotgage[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*age*treat,data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[24]
       #    annowteedis[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*age*treat,data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[9]
       #    annotgeedis[i+16] <- as.data.frame(pairs(emmeans(aov(slope~gen*age*treat,data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[6]
       # }
        
        # stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[2:76]), corr = FALSE,
        #                                                              cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                              symbols = c("****", "***", "**", "*", "")))),
        #                             page=c('Age',as.character(symnum(as.numeric(annoage[2:76]), corr = FALSE,
        #                                                                  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                  symbols = c("****", "***", "**", "*", "")))),
        #                             pgenage=c('Gen:Age',as.character(symnum(as.numeric(annogenage[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pagetreat=c('Age:Treat',as.character(symnum(as.numeric(annoagetreat[2:76]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))),
        #                             pgenagetreat=c('Gen:Age:Treat',as.character(symnum(as.numeric(annogenagetreat[2:76]), corr = FALSE,
        #                                                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                         symbols = c("****", "***", "**", "*", "")))),
        #                             pwtage=c('WT-age',as.character(symnum(as.numeric(annowtage[2:76]), corr = FALSE,
        #                                                                 cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                 symbols = c("****", "***", "**", "*", "")))),
        #                             ptgage=c('TG-age',as.character(symnum(as.numeric(annotgage[2:76]), corr = FALSE,
        #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                     symbols = c("****", "***", "**", "*", "")))),
        #                             pwteedis=c('WT-EEdis',as.character(symnum(as.numeric(annowteedis[2:76]), corr = FALSE,
        #                                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                   symbols = c("****", "***", "**", "*", "")))),
        #                             ptgeedis=c('TG-EEdis',as.character(symnum(as.numeric(annotgeedis[2:76]), corr = FALSE,
        #                                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #                                                                   symbols = c("****", "***", "**", "*", "")))))
      }
      df <- subset(df,df$time!=0)
      # stat.pairwise <- stat.pairwise[c(1:15,17:76),] 
      
      
      if(input$variablecomp4 != "Age effect"){
        p <- ggline(df, x = "time", y = "slope", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=.5,
                    plot_type="p",
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    palette = "Set1",
                    ggtheme = theme_bw(),
                    # title='Input Output Relation',
                    xlab='Time (min.)',
                    ylab='Normalized fEPSP slope',
                    # xlim=c(0,100),
                    ylim=c(75,245),
                    legend='right')  +
          # stat_pvalue_manual(data=stat.pairwise,
          #                    # x=seq(-20,60,1),
          #                    label = "p.signif",
          #                    label.y=rep(210,40),
          #                    xmin=NULL,
          #                    xmax=NULL,
          #                    # hide.ns = T,
          #                    remove.bracket=T,
          #                    size=8
          # ) +
          # stat_compare_means(aes(group = genotype),
        #                    method="anova",
        #                    label = "p.signif",
        #                    label.y=rep(210,40),
        #                    hide.ns = T,
        #                    size=8) +
        theme(axis.text=element_text(size=20),
              axis.title = element_text(size=25),
              title=element_text(size=25),
              legend.text = element_text(size=25))
        # p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
        p #+ 
          # geom_signif(annotations=stat.pairwise$pgentreat,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 245,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pgen,textsize=4,aes(group = genotype),
          #           xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #           y_position = 240,size=0,color='orange3') +
          # geom_signif(annotations=stat.pairwise$ptreat,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 235,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwttg,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 225,size=0,color='#377eb8') +
          # geom_signif(annotations=stat.pairwise$pwtwtee,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 220,size=0,color='#4daf4a') +
          # geom_signif(annotations=stat.pairwise$pwttgee,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 210,size=0,color='#f781bf') +
          # geom_signif(annotations=stat.pairwise$ptgtgee,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 215,size=0,color='#984ea3')  
      }else{
        p <- ggline(df, x = "time", y = "slope", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=.5,
                    plot_type="p",
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    palette = "Set2",
                    ggtheme = theme_bw(),
                    # title='Input Output Relation',
                    xlab='Time (min.)',
                    ylab='Normalized fEPSP slope',
                    # xlim=c(0,100),
                    ylim=c(75,245),
                    legend='right')  +
          theme(axis.text=element_text(size=20),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25))
        p #+ 
          # geom_signif(annotations=stat.pairwise$pgenagetreat,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 240,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pgenage,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 235,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pagetreat,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 230,size=0,color='orange3') +
          # geom_signif(annotations=stat.pairwise$page,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 225,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwteedis,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 222,size=0,color='#66c2a5') +
          # geom_signif(annotations=stat.pairwise$ptgeedis,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 218,size=0,color='#8da0cb') +
          # geom_signif(annotations=stat.pairwise$pwtage,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 210,size=0,color='#66c2a5') +
          # geom_signif(annotations=stat.pairwise$ptgage,textsize=4,aes(group = genotype),
          #             xmin=c(-15:-1,1:60),xmax=c(-15:-1,1:60),tip_length=0,
          #             y_position = 215,size=0,color='#8da0cb')
      }
        
      # print(
      #   ggplotly(p)
      # )
    })
    
    # output$slopeltp <-  renderPlotly({
    output$slopeltp <-  renderPlot({
      print(plotslopeltp())
    })
  
    
    output$downloadPlot3 <- downloadHandler(
      filename = function(){paste('LTP_',input$variablecomp4, '.pdf', sep = '')},
      
      content = function(file){
        ggsave(file,width=12,height=6,plotslopeltp(),useDingbats=F)
      },
      
      contentType = "application/pdf"
    )
    
    # ################
    plotsholl <- reactive({
      
      df = read.csv('sholl_2mo.csv', header = FALSE)
      df2 = data.frame(t(df)) #only for sholl
  
      df3 <- arrange(df2)
      colnames(df3) <- c(rep('WTNE',18),
                         rep('TGNE',19),
                         rep('WTEE',19),
                         rep('TGEE',15))
  
      df3 <- melt(t(df3))
      df3$Var2 <- floor(as.numeric(df3$Var2)*40/2)
      df3_1 <- df3
  
      df = read.csv('sholl_6mo.csv', header = FALSE)
      df2 = data.frame(t(df)) #only for sholl
  
      df3 <- arrange(df2)
      colnames(df3) <- c(rep('WTNE_NE',18),
                         rep('TGNE_NE',18),
                         rep('WTEE_NE',18),
                         rep('TGEE_NE',18))
  
      df3 <- melt(t(df3))
      df3$Var2 <- floor(as.numeric(df3$Var2)*40/2)
  
      df3 <- rbind(df3_1,df3)
      
      blocking <- read.csv('blocking.csv',header=F)
      blocking2mo <- blocking[1:71,]
      blocking6mo <- blocking[72:143,]
      
      df <-  df3
      df$blocking <- as.factor(c(rep(blocking2mo,22),rep(blocking6mo,17)))
      df <- rename(df,genotype=Var1)
      df <- rename(df,radius=Var2)
      df$radius <- as.numeric(df$radius)
      
      if(input$variableshollgroup == "EE 2mo"){
        df <- subset(df,df$genotype=='WTNE'|df$genotype=='WTEE'|df$genotype=='TGNE'|df$genotype=='TGEE')
      }
      else if(input$variableshollgroup == "EE 6mo"){
        df <- subset(df,df$genotype=='WTNE_NE'|df$genotype=='WTEE_NE'|df$genotype=='TGNE_NE'|df$genotype=='TGEE_NE')
        levels(df$genotype)[levels(df$genotype)=="WTNE_NE"] <- "WTNE"
        levels(df$genotype)[levels(df$genotype)=="TGNE_NE"] <- "TGNE"
        levels(df$genotype)[levels(df$genotype)=="WTEE_NE"] <- "WTEE"
        levels(df$genotype)[levels(df$genotype)=="TGEE_NE"] <- "TGEE"
      }
      else{
        # df <- subset(df,df$genotype=='WTNE'|df$genotype=='WTNE_NE'|df$genotype=='TGNE'|df$genotype=='TGNE_NE')
        
        df$gen <- factor('WT',levels=c('WT','TG'))
        df$age <- '6mo'
        df$treat <- NULL
        df$genotype <-  factor(df$genotype,
                               levels=c('WTNE',
                                        'TGNE',
                                        'WTEE',
                                        'TGEE',
                                        'WTNE_NE',
                                        'WTEE_NE',
                                        'TGNE_NE',
                                        'TGEE_NE'))
        df$gen[df$genotype=="WTNE"] <- 'WT'
        df$gen[df$genotype=="WTEE"] <- 'WT'
        df$gen[df$genotype=="TGNE"] <- 'TG'
        df$gen[df$genotype=="TGEE"] <- 'TG'
        df$treat[df$genotype=="WTNE"] <- 'NE'
        df$treat[df$genotype=="WTEE"] <- 'EE'
        df$treat[df$genotype=="TGNE"] <- 'NE'
        df$treat[df$genotype=="TGEE"] <- 'EE'
        df$gen[df$genotype=="WTNE_NE"] <- 'WT'
        df$gen[df$genotype=="WTEE_NE"] <- 'WT'
        df$gen[df$genotype=="TGNE_NE"] <- 'TG'
        df$gen[df$genotype=="TGEE_NE"] <- 'TG'
        df$treat[df$genotype=="WTNE_NE"] <- 'NE'
        df$treat[df$genotype=="WTEE_NE"] <- 'EE'
        df$treat[df$genotype=="TGNE_NE"] <- 'NE'
        df$treat[df$genotype=="TGEE_NE"] <- 'EE'
        df$age[df$genotype=="WTNE" | df$genotype=="TGNE"| df$genotype=="WTEE" | df$genotype=="TGEE"] <- '2mo'
        
        levels(df$genotype)[levels(df$genotype)=="WTNE"] <- 'WTNE_2mo'
        levels(df$genotype)[levels(df$genotype)=="TGNE"] <- 'TGNE_2mo'
        levels(df$genotype)[levels(df$genotype)=="WTEE"] <- 'WTEE'
        levels(df$genotype)[levels(df$genotype)=="TGEE"] <- 'TGEE'
        levels(df$genotype)[levels(df$genotype)=="WTNE_NE"] <- 'WTNE_6mo'
        levels(df$genotype)[levels(df$genotype)=="TGNE_NE"] <- 'TGNE_6mo'
        levels(df$genotype)[levels(df$genotype)=="WTEE_NE"] <- 'WTEEdis'
        levels(df$genotype)[levels(df$genotype)=="TGEE_NE"] <- 'TGEEdis'
        
      }
      if(input$variableshollgroup != "Age effect"){
        df$genotype <- factor(df$genotype,
                              levels=c('WTNE',
                                       'TGNE',
                                       'WTEE',
                                       'TGEE'))
      }else{
        df$genotype <- factor(df$genotype,levels=c('WTNE_2mo',
                                                     'WTNE_6mo',
                                                     'TGNE_2mo',
                                                     'TGNE_6mo',
                                                     'WTEE',
                                                     'WTEEdis',
                                                     'TGEE',
                                                     'TGEEdis'))
        df$alphav <- NULL
        df$alphav[df$genotype %in% c('WTNE_2mo','TGNE_2mo','WTEE','TGEE')] <- 1
        df$alphav[df$genotype %in% c('WTNE_6mo','TGNE_6mo','WTEEdis','TGEEdis')] <- 0.6
      }
      
      # if(input$variableshollgroup != "Age effect"){
      #   df$gen <- sapply(df$genotype,substring,first=1,last=2)
      #   df$treat <- sapply(df$genotype,substring,first=3,last=4)
      #   
      #   annogen <- NULL
      #   annotreat <- NULL
      #   annogentreat <- NULL
      #   annowttg <- NULL
      #   annowtwtee <- NULL
      #   annowttgee <- NULL
      #   annotgtgee <- NULL
      #   for(i in seq(20,440,20)){
      #     dfi <- subset(df,df$radius==i)
      #     tryCatch(print(summary(mixed(value~gen*treat+(1|blocking),data=dfi))), error=function(e){1})
      #     # tryCatch(print(summary(aov(value~gen*treat,data=dfi))), error=function(e){})
      #     print(i)
      #     groupannop <- tryCatch(as.data.frame(summary(mixed(value~gen*treat+(1|blocking),data=dfi))[[10]])[2:4,5],error=function(e){c(1,1,1)})
      #     # groupannop <- tryCatch(as.data.frame(summary(aov(value~gen*treat,data=dfi))[[10]])[2:4,5],error=function(e){c(1,1,1)})
      #     
      #     annogen[i/20] <- groupannop[1]
      #     annotreat[i/20] <- groupannop[2]
      #     annogentreat[i/20] <- groupannop[3]
      #     if(i==140|i==160) dfi2 <-dfi
      #     if(i==140|i==160)tryCatch(print(pairs(emmeans(mixed(value~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm")),error=function(e){})
      #     annowttg[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[6],error=function(e){1})
      #     annowtwtee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[5],error=function(e){1})
      #     annowttgee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[3],error=function(e){1})
      #     annotgtgee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*treat+(1|blocking),data=dfi), ~gen*treat), adjust = "holm"))$p.value[2],error=function(e){1})
      #     
      #     # annowttg[i/20] <- tryCatch(as.data.frame(pairs(emmeans(aov(value~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[6],error=function(e){1})
      #     # annowtwtee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(aov(value~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[5],error=function(e){1})
      #     # annowttgee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(aov(value~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[3],error=function(e){1})
      #     # annotgtgee[i/20] <- tryCatch(as.data.frame(pairs(emmeans(aov(value~gen*treat,data=dfi), ~gen*treat), adjust = "holm"))$p.value[2],error=function(e){1})
      #   }
      #   
      #   annogen[is.na(annogen[1:22])] <- 1
      #   annotreat[is.na(annotreat[1:22])] <- 1
      #   annogentreat[is.na(annogentreat[1:22])] <- 1
      #   annowttg[is.na(annowttg[1:22])] <- 1
      #   annowtwtee[is.na(annowtwtee[1:22])] <- 1
      #   annowttgee[is.na(annowttgee[1:22])] <- 1
      #   annotgtgee[is.na(annotgtgee[1:22])] <- 1
      #   
      #   stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[1:22]), corr = FALSE,
      #                                                                cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                symbols = c("****", "***", "**", "*", "")))),
      #                               ptreat=c('Treat',as.character(symnum(as.numeric(annotreat[1:22]), corr = FALSE,
      #                                                                    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                    symbols = c("****", "***", "**", "*", "")))),
      #                               pgentreat=c('Gen:Treat',as.character(symnum(as.numeric(annogentreat[1:22]), corr = FALSE,
      #                                                                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                           symbols = c("****", "***", "**", "*", "")))),
      #                               pwttg=c('WT-TG',as.character(symnum(as.numeric(annowttg[1:22]), corr = FALSE,
      #                                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                   symbols = c("****", "***", "**", "*", "")))),
      #                               pwtwtee=c('WT-WTEE',as.character(symnum(as.numeric(annowtwtee[1:22]), corr = FALSE,
      #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                       symbols = c("****", "***", "**", "*", "")))),
      #                               pwttgee=c('WT-TGEE',as.character(symnum(as.numeric(annowttgee[1:22]), corr = FALSE,
      #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                       symbols = c("****", "***", "**", "*", "")))),
      #                               ptgtgee=c('TG-TGEE',as.character(symnum(as.numeric(annotgtgee[1:22]), corr = FALSE,
      #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                       symbols = c("****", "***", "**", "*", "")))))
      # }else{
      #   # df$gen <- sapply(df$genotype,substring,first=1,last=2)
      #   # df$age <- sapply(df$genotype,substring,first=3,last=5)
      #   
      #   annogen <- NULL
      #   annoage <- NULL
      #   annogenage <- NULL
      #   annoagetreat <- NULL
      #   annodisc <- NULL
      #   annowtage <- NULL
      #   annotgage <- NULL
      #   annowteedis <- NULL
      #   annotgeedis <- NULL
      #   for(i in seq(20,440,20)){
      #     dfi <- subset(df,df$radius==i)
      #     tryCatch(print(summary(mixed(value~gen*age*treat+(1|blocking),data=dfi))), error=function(e){1})
      #     print(i)
      #     groupannop <- tryCatch(as.data.frame(summary(mixed(value~gen*age*treat+(1|blocking),data=dfi))[[10]])[2:7,5],error=function(e){c(1,1,1)})
      #     # groupannop <- as.data.frame(summary(mixed(value~gen*age+(1|blocking),data=dfi))[[10]])[2:4,5]
      #     annogen[i/20] <- groupannop[1]
      #     agenep <- tryCatch(as.data.frame(joint_tests(mixed(value~gen*age*treat+(1|blocking),data=dfi), by = "treat"))$p.value[4],error=function(e){1})
      #     annoage[i/20] <- agenep  #groupannop[2]
      #     annogenage[i/20] <- groupannop[3]
      #     annoagetreat[i/20] <- groupannop[6]
      #     # print(as.data.frame(joint_tests(mixed(values~gen*age*treat+(1|blocking),data=dfi), by = "treat")))
      #     discannop <- tryCatch(as.data.frame(joint_tests(mixed(value~gen*age*treat+(1|blocking),data=dfi), by = "treat"))$p.value[3],error=function(e){1})
      #     annodisc[i/20] <- discannop
      #     annowtage[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[27],error=function(e){1})
      #     tryCatch(print(pairs(emmeans(mixed(value~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm")),error=function(e){})
      #     annotgage[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[24],error=function(e){1})
      #     annowteedis[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[9],error=function(e){1})
      #     annotgeedis[i/20] <- tryCatch(as.data.frame(pairs(emmeans(mixed(value~gen*age*treat+(1|blocking),data=dfi), ~gen*age*treat), adjust = "holm"))$p.value[2],error=function(e){1})
      #   }
      #   
      #   annogen[is.na(annogen[1:22])] <-1
      #   annoage[is.na(annoage[1:22])] <-1
      #   annogenage[is.na(annogenage[1:22])] <-1
      #   annoagetreat[is.na(annoagetreat[1:22])] <-1
      #   annodisc[is.na(annodisc[1:22])] <-1
      #   annowtage[is.na(annowtage[1:22])] <-1
      #   annotgage[is.na(annotgage[1:22])] <-1
      #   annowteedis[is.na(annowteedis[1:22])] <-1
      #   annotgeedis[is.na(annotgeedis[1:22])] <-1
      #   
      #   stat.pairwise <- data.frame(pgen=c('Gen',as.character(symnum(as.numeric(annogen[1:22]), corr = FALSE,
      #                                                                cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                symbols = c("****", "***", "**", "*", "")))),
      #                               page=c('Age NE',as.character(symnum(as.numeric(annoage[1:22]), corr = FALSE,
      #                                                                cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                symbols = c("****", "***", "**", "*", "")))),
      #                               pgenage=c('Gen:Age',as.character(symnum(as.numeric(annogenage[1:22]), corr = FALSE,
      #                                                                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                       symbols = c("****", "***", "**", "*", "")))),
      #                               # pagetreat=c('Age:Treat',as.character(symnum(as.numeric(annoagetreat[1:22]), corr = FALSE,
      #                               #                                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                               #                                         symbols = c("****", "***", "**", "*", "")))),
      #                               pagetreat=c('Disc',as.character(symnum(as.numeric(annodisc[1:22]), corr = FALSE,
      #                                                                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                           symbols = c("****", "***", "**", "*", "")))),
      #                               pwtage=c('WT-age',as.character(symnum(as.numeric(annowtage[1:22]), corr = FALSE,
      #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                     symbols = c("****", "***", "**", "*", "")))),
      #                               ptgage=c('TG-age',as.character(symnum(as.numeric(annotgage[1:22]), corr = FALSE,
      #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                     symbols = c("****", "***", "**", "*", "")))),
      #                               pwteedis=c('WT-EEdis',as.character(symnum(as.numeric(annowteedis[1:22]), corr = FALSE,
      #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                     symbols = c("****", "***", "**", "*", "")))),
      #                               ptgeedis=c('TG-EEdis',as.character(symnum(as.numeric(annotgeedis[1:22]), corr = FALSE,
      #                                                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                                                     symbols = c("****", "***", "**", "*", "")))))
      # }
      # df <- subset(df,df$time!=0)
      # stat.pairwise <- stat.pairwise[1:11,]
      
      
      if(input$variableshollgroup != "Age effect"){
        p <- ggline(df, x = "radius", y = "value", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=1.5,
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    palette = "Set1",
                    ggtheme = theme_pubr(),
                    # title='Input Output Relation',
                    xlab='Radius [um]',
                    ylab='Number of intersections',
                    xlim=c(-1,22),
                    ylim=c(0,18),
                    legend='right') +
          theme(axis.text=element_text(size=15),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25)) 
        # p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
        
        p #+ 
          # geom_signif(annotations=stat.pairwise$pgentreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 17,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pgen,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 16,size=0,color='orange3') +
          # geom_signif(annotations=stat.pairwise$ptreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 15,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwttg,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 14,size=0,color='#377eb8') +
          # geom_signif(annotations=stat.pairwise$pwtwtee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 13,size=0,color='#4daf4a') +
          # geom_signif(annotations=stat.pairwise$pwttgee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 12,size=0,color='#f781bf') +
          # geom_signif(annotations=stat.pairwise$ptgtgee,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 11,size=0,color='#984ea3')  
      }else{
        palettev <- c("#E41A1C","#E41A1C","#377EB8","#377EB8","#4DAF4A","#4DAF4A","#984EA3","#984EA3")
        # df <- subset(df,df$genotype %in% c('WTEE','WTEEdis','TGEE','TGEEdis'))
        p <- ggline(df, x = "radius", y = "value", add = c("mean_se"),#,"point"),
                    color = "genotype", 
                    size=1.5,
                    add.params = list(alpha=0.5, size=0.7),
                    shape = "genotype",
                    linetype = "genotype",
                    # palette = "8-class Paired",
                    palette = palettev,
                    alpha = "alphav",
                    ggtheme = theme_pubr(),
                    # title='Input Output Relation',
                    xlab='Radius [um]',
                    ylab='Number of intersections',
                    xlim=c(-1,22),
                    ylim=c(0,16.5),
                    legend='right') +
          theme(axis.text=element_text(size=15),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25)) 
        # p$layers <- c(stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),p$layers)
        p #+ 
          # geom_signif(annotations=stat.pairwise$page,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 16,size=0,color='brown') +
          # geom_signif(annotations=stat.pairwise$pagetreat,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 15,size=0,color='orange3') +
          # # geom_signif(annotations=stat.pairwise$page,textsize=4,aes(group = genotype),
          # #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          # #             y_position = 11,size=0,color='orange') +
          # geom_signif(annotations=stat.pairwise$pwtage,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 11,size=0,color='#fb9a99') +
          # geom_signif(annotations=stat.pairwise$ptgage,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 12,size=0,color='#b2df8a') +
          # geom_signif(annotations=stat.pairwise$pwteedis,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 13,size=0,color='#a6cee3') +
          # geom_signif(annotations=stat.pairwise$ptgeedis,textsize=4,aes(group = genotype),
          #             xmin=seq(0,22),xmax=seq(0,22),tip_length=0,
          #             y_position = 14,size=0,color='#e78ac3')
      }
      
      # print(
      #   ggplotly(p)
      # )
    })
    
    output$sholl <-  renderPlot({
      # print(plotsholl())
      plotsholl()
    })
    
    
    output$downloadPlotsholl <- downloadHandler(
      filename = function(){paste('Sholl_',input$variableshollgroup, '.pdf', sep = '')},
      
      content = function(file){
        ggsave(file,width=12,height=8,plotsholl(),useDingbats=F)
      },
      
      contentType = "application/pdf"
    )
    # output$sholl <- renderPlot({
    #   df = read.csv('sholl_2mo.csv', header = FALSE)
    #   df2 = data.frame(t(df)) #only for sholl
    #   
    #   df3 <- arrange(df2)
    #   colnames(df3) <- c(rep('WT',18),
    #                      rep('TgDyrk1A',19),
    #                      rep('WT-EE',19),
    #                      rep('TgDyrk1A-EE',15))
    #   
    #   df3 <- melt(t(df3))
    #   df3$Var2 <- floor(df3$Var2*40/2)
    #   df3_1 <- df3
    #   
    #   df = read.csv('sholl_6mo.csv', header = FALSE)
    #   df2 = data.frame(t(df)) #only for sholl
    #   
    #   df3 <- arrange(df2)
    #   colnames(df3) <- c(rep('WT-NE_NE',18),
    #                      rep('TgDyrk1A-NE_NE',18),
    #                      rep('WT-EE_NE',18),
    #                      rep('TgDyrk1A-EE_NE',18))
    #   
    #   df3 <- melt(t(df3))
    #   df3$Var2 <- floor(df3$Var2*40/2)
    #   
    #   df3 <- rbind(df3_1,df3)
    #   
    #   df3 <- filter(df3, Var1 %in% input$variablegroups)
    #   
    #   df3 <- plyr::rename(df3, c("Var1" = "Genotype"))
    #   
    #   ggline(df3, x = "Var2", y = "value", add = "mean_se", 
    #          size=1.5,
    #          color = "Genotype", 
    #          shape = "Genotype",
    #          linetype = "Genotype",
    #          palette = "Set1",
    #          # title='Sholl EE 2mo',
    #          xlab='Distance to soma [micrometers]',
    #          ylab='Number of intersections',
    #          legend='right')+
    #     stat_compare_means(aes(group = Genotype), 
    #                        method="anova", 
    #                        label = "p.signif", 
    #                        label.y=12,
    #                        size=8,
    #                        hide.ns=T) +
    #     theme(axis.text=element_text(size=20),
    #           axis.title = element_text(size=25),
    #           title=element_text(size=25),
    #           legend.text = element_text(size=25)) +
    #     scale_x_discrete(breaks=seq(20,440,40))
    #     
    #   
    # })
    #############################
    # AP layer width
    # output$width <-  renderPlot({
    plotwidth <-  reactive({
        
      if(input$variablecomp=="EE 2mo"){
        data = read.csv('Widths_CA1_2mo.csv')
        data$genotype <- factor(data$genotype,levels=c('WT','TG','WT-EE','TG-EE'))
        # levels(data$genotype)[levels(data$genotype)=="WT"] <- "Wild-Type"
        levels(data$genotype)[levels(data$genotype)=="TG"] <- "TgDyrk1A"
        levels(data$genotype)[levels(data$genotype)=="WT-EE"] <- "WT_EE"
        levels(data$genotype)[levels(data$genotype)=="TG-EE"] <- "TgDyrk1A_EE"
        # data$med_cent <- factor(data$med_cent,levels=c('medial','central'))
        # data = subset(data,data$med_cent==2)
        data = subset(data,data$slanted_jl==0)
        # data = subset(data,data$slanted_jl==0 & data$bregma>)
      }
      else if(input$variablecomp=="EE 6mo"){
        data = read.csv('Widths_CA1_6mo.csv')
        data$genotype <- factor(data$genotype,levels=c('WT-NE-NE','TG-NE-NE','WT-EE-NE','TG-EE-NE'))
        levels(data$genotype)[levels(data$genotype)=="WT-NE-NE"] <- "WT"
        levels(data$genotype)[levels(data$genotype)=="TG-NE-NE"] <- "TgDyrk1A"
        levels(data$genotype)[levels(data$genotype)=="WT-EE-NE"] <- "WT_EE"
        levels(data$genotype)[levels(data$genotype)=="TG-EE-NE"] <- "TgDyrk1A_EE"
        # data$med_cent <- factor(data$med_cent,levels=c('medial','central'))
        # data = subset(data,data$med_cent==2)
        data = subset(data,data$slanted_jl==0)
      }
      else if(input$variablecomp=="Age effect"){
        data1 = read.csv('Widths_CA1_2mo.csv')
        data1$gen <- NULL
        data1$age <- NULL
        data1$treat <- NULL
        data1$gen[data1$genotype=="WT"] <- 'WT'
        data1$gen[data1$genotype=="WT-EE"] <- 'WT'
        data1$gen[data1$genotype=="TG"] <- 'TG'
        data1$gen[data1$genotype=="TG-EE"] <- 'TG'
        data1$treat[data1$genotype=="WT"] <- 'NE'
        data1$treat[data1$genotype=="WT-EE"] <- 'EE'
        data1$treat[data1$genotype=="TG"] <- 'NE'
        data1$treat[data1$genotype=="TG-EE"] <- 'EE'
        data1$age <- '2mo'
        
        data1$genotype <- as.factor(data1$genotype)
        # data1 <- subset(data1,data1$genotype=='WTNE' | data1$genotype=='TGNE') 
        levels(data1$genotype)[levels(data1$genotype)=="WT"] <- 'WTNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="TG"] <- 'TGNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="WT-EE"] <- 'WTEE'
        levels(data1$genotype)[levels(data1$genotype)=="TG-EE"] <- 'TGEE'
        
        data2 = read.csv('Widths_CA1_6mo.csv')
        data2$gen <- NULL
        data2$age <- NULL
        data2$treat <- NULL
        data2$gen[data2$genotype=="WT-NE-NE"] <- 'WT'
        data2$gen[data2$genotype=="WT-EE-NE"] <- 'WT'
        data2$gen[data2$genotype=="TG-NE-NE"] <- 'TG'
        data2$gen[data2$genotype=="TG-EE-NE"] <- 'TG'
        data2$treat[data2$genotype=="WT-NE-NE"] <- 'NE'
        data2$treat[data2$genotype=="WT-EE-NE"] <- 'EE'
        data2$treat[data2$genotype=="TG-NE-NE"] <- 'NE'
        data2$treat[data2$genotype=="TG-EE-NE"] <- 'EE'
        data2$age <- '6mo'
        data2$genotype <- as.factor(data2$genotype)
        # data2 <- subset(data2,data2$genotype=='WTNE' | data2$genotype=='TGNE') 
        levels(data2$genotype)[levels(data2$genotype)=="WT-NE-NE"] <- 'WTNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="TG-NE-NE"] <- 'TGNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="WT-EE-NE"] <- 'WTEEdis'
        levels(data2$genotype)[levels(data2$genotype)=="TG-EE-NE"] <- 'TGEEdis'
        
        data1$bregma <- NULL
        data <- rbind(data1,data2)
        
        # data$genotype <- factor(data$genotype,levels=c('Wild-Type_2mo','Wild-Type_6mo','TgDyrk1A_2mo','TgDyrk1A_6mo'))
        data$genotype <- factor(data$genotype,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
        
        data$alphav <- NULL
        data$alphav[data$genotype %in% c('WTNE_2mo','TGNE_2mo','WTEE','TGEE')] <- 1
        data$alphav[data$genotype %in% c('WTNE_6mo','TGNE_6mo','WTEEdis','TGEEdis')] <- 0.6
        levels(data$alphav) <- seq(0,1,0.1)
      }
      
      data$metric <- data$str_rad + data$str_lac
      data$blocking <- data$Animal
      data$group <- data$genotype
      df <- data
      
      bvdplot <- function(df, metric_name) { 
        if(input$variablecomp=="Age effect"){
          df$group <- factor(df$group,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
          # print(unique(df$group))
          
          df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking),gen=df$gen,age=df$age,treat=df$treat)
          print(summary(mixed(values~gen*age*treat+(1|blocking),data=df_st)))
          print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "age")))
          print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "gen")))
          print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "treat")))
          print(as.data.frame(pairs(emmeans(mixed(values~gen*age*treat+(1|blocking),data=df_st), ~gen*age*treat), adjust = "holm")))
          anno <- as.data.frame(pairs(emmeans(mixed(values~gen*age*treat+(1|blocking),data=df_st), ~gen*age*treat), adjust = "holm"))$p.value
          
          # Age
          # print("Age")
          # print(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]))
          # print(summary(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),])))
          # print(pairs(emmeans(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]), ~gen*age), adjust = "holm"))
          # # Discontinuation
          # print("Discontinuation")
          # print(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]))
          # print(summary(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),])))
          # print(pairs(emmeans(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]), ~gen*age), adjust = "holm"))
          
        }
        else{
          df$group <- factor(df$group,levels=c('WT','TgDyrk1A','WT_EE','TgDyrk1A_EE'))
        
          # Kruskal test to take then format
          stat.all <- compare_means(
            metric ~ group, data = df,
            method = "kruskal.test"
          )
          
          # Pairwise wilcox-test between groups to take then format
          stat.pairwise <- compare_means(
            metric ~ group, data = df,
            method = "wilcox.test"
          ) %>%
            mutate(y.position = max(df$metric)*seq(1.05,1.55,length.out=6))
          
          stat.pairwise$y.position[6] <- stat.pairwise$y.position[4]
          stat.pairwise$y.position[5] <- stat.pairwise$y.position[4]
          print('stat.pairwise')
          print(stat.pairwise)
          stat.pairwise <- stat.pairwise[c(1:3,5),]
          if(input$variablecomp=="Age effect"){
            stat.pairwise <- stat.pairwise[c(1:3,6),]
          }
          
          # Linear mixed effects model
          df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
          # mod <- lmer(values~Genotypes+(1|blocking),data=df_st)
          # coefs <- data.frame(coef(summary(mod)))
          # get p valueues
          # df.KR <- get_ddf_Lb(mod, fixef(mod))
          # p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
          # coefs$p.value <- p.KR
          # LMM using mixed function in afex package
          m1 <- mixed(values~Genotypes+(1|blocking),data=df_st)
          print(summary(m1))
          full <- nice(m1)
          # print(full)
          #print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
          (emm1 <- emmeans(m1, "Genotypes"))
          prs <- as.data.frame(pairs(emm1, adjust = "holm"))
          print(prs)
          
          #coefs
          shap_wtne <- shapiro.test(df$metric[df$group==unique(df$group)[1]])
          shap_tgne <- shapiro.test(df$metric[df$group==unique(df$group)[2]])
          shap_wtee <- shapiro.test(df$metric[df$group==unique(df$group)[3]])
          shap_tgee <- shapiro.test(df$metric[df$group==unique(df$group)[4]])
          annolmm <- as.data.frame(full)
          annolmm[1,3] <- strsplit(annolmm[1,3],' ')[1]
          annolmm[1,4] <- as.numeric(annolmm[1,4])
          annolmm <- rename(annolmm,'contrast'='Effect')
          annolmm <- rename(annolmm,'F/t.ratio'='F')
          annoprs <- prs[c(1:3,5),c(1,4,5,6)]
          if(input$variablecomp=="Age effect"){
            annoprs <- prs[c(1:3,6),c(1,4,5,6)]
          }
          annoprs$df <- round(annoprs$df, digits=2)
          annoprs$t.ratio <- round(annoprs$t.ratio, digits=2)
          annoprs$p.value <- formatC(annoprs$p.value,digits=2,format = 'g')
          
          annoprs <- rename(annoprs,'F/t.ratio'='t.ratio')
          print(annolmm)
          print(annoprs)
          anno <- rbind(annolmm,annoprs)
          pshap <- c(round(shap_wtne$p.value, digits = 2),
                     round(shap_tgne$p.value, digits = 2),
                     round(shap_wtee$p.value, digits = 2),
                     round(shap_tgee$p.value, digits = 2),
                     '-')
          metricnm <- rep(metric_name,5)
          group <- c('WT',
                     'TgDyrk1A',
                     'WT_EE',
                     'TgDyrk1A_EE',
                     '-'
          )
          
          means <- c(
            round(mean(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
            round(mean(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
            round(mean(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
            round(mean(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
            '-'
          )
          sds <- c(
            round(sd(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
            round(sd(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
            round(sd(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
            round(sd(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
            '-'
          )
          anno <- cbind(metricnm,group,means,sds,pshap,anno)
          
          print(anno)
          print(xtable(anno, caption = "Statistics for..."), include.rownames = FALSE)
          
          # modify test data frame with linear mixed effects model p values
          stat.all$p <- as.numeric(anno$p.value[1])
          stat.all$p.adj <- anno$p.value[1]
          stat.all$p.format <- substr(as.character(anno$p.value[1]),1,6)
          stat.all$p.signif <- as.character(symnum(as.numeric(stat.all$p), corr = FALSE,
                                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                   symbols = c("****", "***", "**", "*", "ns")))
          stat.all$method <- 'Linear.Mixed-Effects'
          
          stat.pairwise$p <- c(anno$p.value[2],
                               anno$p.value[3],
                               anno$p.value[4],
                               anno$p.value[5]
          )
          
          stat.pairwise$p.adj <- c(anno$p.value[2],
                                   anno$p.value[3],
                                   anno$p.value[4],
                                   anno$p.value[5]
          )
          
          stat.pairwise$p.format <- c(substr(as.character(anno$p.value[2]),1,6),
                                      substr(as.character(anno$p.value[3]),1,6),
                                      substr(as.character(anno$p.value[4]),1,6),
                                      substr(as.character(anno$p.value[5]),1,6)
          )
          
          stat.pairwise$p.signif <- as.character(symnum(as.numeric(stat.pairwise$p), corr = FALSE,
                                                        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                        symbols = c("****", "***", "**", "*", "ns")))
          
          stat.pairwise$method <- rep('Linear.Mixed-Effects',4)
          
          anno <- as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), ~Genotypes), adjust = "holm"))$p.value[c(1,6)]
        }
        
        if(input$variablecomp=="Age effect"){
          # print(df$group)
          palettev <- c("#E41A1C","#E41A1C","#377EB8","#377EB8","#4DAF4A","#4DAF4A","#984EA3","#984EA3")
          # alphav <- c(1,0.4,1,0.4,1,0.4,1,0.4)
          # print(df$alphav)
          ggviolin(df, x="group", y="metric", fill = "group",width=0.5, 
                   trim=TRUE,
                   # palette="8-class Paired",
                   palette = palettev,
                   alpha = "alphav",
                   add = c("boxplot","dotplot"),
                   add.params = list(size=0.5,alpha = "alphav"),
                   ylab = metric_name) #+
            # geom_hline(yintercept=median(df$metric[df$group=="WTNE_2mo"]),linetype=2) +
            # stat_pvalue_manual(
            #   data = stat.pairwise, label = "p.signif",
            #   # xmin = "group1", xmax = "group2",
            #   y.position = "y.position",
            #   label.size = 6#,
            #   # comparisons=list(c(1,2),
            #   #                  c(3,4)
            #   # )
            # ) +
            # geom_signif(annotation=formatC(anno[c(27,24,9,2)], digits=2),
            #             aes(group = group), 
            #             map_signif_level = T,
            #             textsize=6,
            #             y_position = rep(max(df$metric)*1.1,4),
            #             comparisons=list(c(1,2),
            #                              c(3,4),
            #                              c(5,6),
            #                              c(7,8)
            #             )
            # )
        }
        else {
          ggviolin(df, x="group", y="metric", fill = "group",width=0.5, 
                 trim=TRUE,
                 palette="Set1",
                 add = c("boxplot","dotplot"),
                 add.params = list(size=0.5),
                 ylab = metric_name) +
          geom_hline(yintercept=median(df$metric[df$group=="WT"]),linetype=2) +
          theme(axis.text=element_text(size=20),
                  axis.title = element_text(size=25),
                  title=element_text(size=25),
                  legend.text = element_text(size=25)) #+
            # stat_pvalue_manual(
            #   data = stat.pairwise, label = "p.signif",
            #   xmin = "group1", xmax = "group2",
            #   y.position = "y.position",
            #   label.size=8
            # )
        }
      }
      
      # df <- filter(df, group %in% input$variablegroups)
      
      # print(
      #   ggplotly(
          if(input$variablesubl == 'Apical CA1'){
            bvdplot(df,'CA1 width')
          }
          else if(input$variablesubl == 'Stratum Radiatum'){
            df$metric <- df$str_rad
            bvdplot(df,'Stratum Radiatum width')
          }
          else if(input$variablesubl == 'Stratum Lacunosum'){
            df$metric <- df$str_lac
            bvdplot(df,'Stratum Lacunosum width')
          }else if(input$variablesubl == 'Stratum Oriens'){
            df$metric <- df$str_lac
            bvdplot(df,'Stratum Oriens width')
          }
      #   )
      # )
    })
    
    output$width <-  renderPlot({
      print(plotwidth())
    })
    
    output$downloadPlot4 <- downloadHandler(
      filename = function(){paste(input$variablesubl,'_',input$variablecomp, '.pdf', sep = '')},
      
      content = function(file){
        ggsave(file,width=12,height=6,plotwidth(),useDingbats=F)
      },
      
      contentType = "application/pdf"
    )
    
    output$tablewidth <- renderTable({
      if(input$variablecomp == 'EE 2mo'){
        data <- read.csv('Widths_CA1_2mo.csv')
        levels(data$genotype)[levels(data$genotype)=="WT"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TG"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="WT-EE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TG-EE"] <- "TGEE"
        data = subset(data,data$slanted_jl==0)
      }
      else if(input$variablecomp == 'EE 6mo'){
        data <- read.csv('Widths_CA1_6mo.csv')
        levels(data$genotype)[levels(data$genotype)=="WT-NE-NE"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TG-NE-NE"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="WT-EE-NE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TG-EE-NE"] <- "TGEE"
        data = subset(data,data$slanted_jl==0)
      }
      else if(input$variablecomp == 'Age effect'){
        data1 <- read.csv('Widths_CA1_2mo.csv')
        data1 = subset(data1,data1$slanted_jl==0)
        data2 <- read.csv('Widths_CA1_6mo.csv')
        data2 = subset(data2,data2$slanted_jl==0)
        data1$bregma <- NULL
        data=rbind(data1,data2)
        data <- subset(data,data$genotype!="WT-EE"&data$genotype!="TG-EE"&data$genotype!="WT-EE-NE"&data$genotype!="TG-EE-NE")
        levels(data$genotype)[levels(data$genotype)=="WT"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TG"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="WT-NE-NE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TG-NE-NE"] <- "TGEE"
      }
      
      if(input$variablesubl == 'Apical CA1'){
        data$both <- data$str_lac+data$str_rad
        subl <- 'both'
      }
      else if(input$variablesubl == 'Stratum Radiatum'){
        subl <- 'str_rad'
      }
      else if(input$variablesubl == 'Stratum Lacunosum'){
        subl <- 'str_lac'
      }
      else if(input$variablesubl == 'Stratum Oriens'){
        subl <- 'str_or'
      }
      # df <- data
      # df$gen <- sapply(df$genotype,substring,first=1,last=2)
      # df$treat <- sapply(df$genotype,substring,first=3,last=4)
      # df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,treat=df$treat)
      # # df$group <- factor(df$group,levels=c('WT','TgDyrk1A','WT_EE','TgDyrk1A_EE'))
      # df$group <- factor(df$group,levels=c('WTNE','TGNE','WTEE','TGEE'))
      # df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      # 
      if(input$variablecomp != 'Age effect'){
        df <- data
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$treat <- sapply(df$genotype,substring,first=3,last=4)
        df <- data.frame(metric=data[,which(names(data)==subl)],group=data$genotype,blocking=data$Animal,gen=df$gen,treat=df$treat)
        df$group <- factor(df$group,levels=c('WTNE','TGNE','WTEE','TGEE'))
        df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*treat+(1|blocking),data=df))
        print(summary(mixed(metric~gen*treat+(1|blocking),data=df)))
        groupanno <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Treatment effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Treatment interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(6,5,3,2)]
        groupannodf <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.den.Df
        groupannoF <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.F
        groupannop <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
      }
      else{
        df <- data
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$age <- sapply(df$genotype,substring,first=3,last=4)
        df <- data.frame(metric=data[,which(names(data)==subl)],group=data$genotype,blocking=data$Animal,gen=df$gen,age=df$age)
        df$group <- factor(df$group,levels=c('WTNE','WTEE','TGNE','TGEE'))
        levels(df$group)[levels(df$group)=="WTNE"] <- "Wild-Type_2mo"
        levels(df$group)[levels(df$group)=="TGNE"] <- "TgDyrk1A_2mo"
        levels(df$group)[levels(df$group)=="WTEE"] <- "Wild-Type_6mo"
        levels(df$group)[levels(df$group)=="TGEE"] <- "TgDyrk1A_6mo"
        df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*age+(1|blocking),data=df))
        groupanno <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Age effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Age interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*age+(1|blocking),data=df), ~gen*age), adjust = "holm"))$p.value[c(5,2)]
        groupannodf <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.den.Df
        groupannoF <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.F
        groupannop <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.Pr..F.
      }
      # anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(6,5,3,2)]
      
      
      
      # Kruskal test to take then format
      stat.all <- compare_means(
        metric ~ group, data = df,
        method = "kruskal.test"
      )
      
      # Pairwise wilcox-test between groups to take then format
      stat.pairwise <- compare_means(
        metric ~ group, data = df,
        method = "wilcox.test"
      ) %>%
        mutate(y.position = max(df$metric)*seq(1.05,1.55,length.out=6))
      
      stat.pairwise$y.position[6] <- stat.pairwise$y.position[4]
      stat.pairwise$y.position[5] <- stat.pairwise$y.position[4]
      stat.pairwise <- stat.pairwise[c(1:3,5),]
      if(input$variablecomp=="Age effect"){
        stat.pairwise <- stat.pairwise[c(1:3,6),]
      }
      
      # Linear mixed effects model
      df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      m1 <- mixed(values~Genotypes+(1|blocking),data=df_st)
      # full <- nice(m1)
      # print(full)
      #print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
      (emm1 <- emmeans(m1, "Genotypes"))
      prs <- as.data.frame(pairs(emm1, adjust = "holm"))
      print(prs)
      
      #coefs
      shap_wtne <- shapiro.test(df$metric[df$group==unique(df$group)[1]])
      shap_tgne <- shapiro.test(df$metric[df$group==unique(df$group)[2]])
      shap_wtee <- shapiro.test(df$metric[df$group==unique(df$group)[3]])
      shap_tgee <- shapiro.test(df$metric[df$group==unique(df$group)[4]])
      # annolmm <- as.data.frame(full)
      # annolmm[1,3] <- strsplit(annolmm[1,3],' ')[1]
      # annolmm[1,4] <- as.numeric(annolmm[1,4])
      # annolmm <- rename(annolmm,'contrast'='Effect')
      # annolmm <- rename(annolmm,'F/t.ratio'='F')
      annoprs0 <- prs[c(1:3,5),c(1,4,5,6)]
      if(input$variablecomp=="Age effect"){
        annoprs0 <- prs[c(1:3,6),c(1,4,5,6)]
      }
      # annoprs$df <- round(annoprs$df, digits=2)
      # annoprs$t.ratio <- round(annoprs$t.ratio, digits=2)
      # annoprs$p.value <- formatC(annoprs$p.value,digits=2,format = 'g')
      annoprs <- data.frame()[1:7,]
      annoprs$df <- c(formatC(groupannodf,digits=2),round(annoprs0$df, digits=2))
      annoprs$t.ratio <- c(formatC(groupannoF,digits=2),round(annoprs0$t.ratio, digits=2))
      annoprs$p.value <- c(formatC(groupannop,digits=2),formatC(annoprs0$p.value,digits=2,format = 'g'))
      
      annoprs <- rename(as.data.frame(annoprs),'F/t.ratio'='t.ratio')
      # print(annolmm)
      # print(annoprs)
      # anno <- rbind(annolmm,annoprs)
      
      pshap <- c(round(shap_wtne$p.value, digits = 2),
                 round(shap_tgne$p.value, digits = 2),
                 round(shap_wtee$p.value, digits = 2),
                 round(shap_tgee$p.value, digits = 2),
                 '',
                 '',
                 '')
      if(input$variablecomp=="Age effect"){
        pshap <- c(round(shap_wtne$p.value, digits = 2),
                   round(shap_tgne$p.value, digits = 2),
                   round(shap_wtee$p.value, digits = 2),
                   round(shap_tgee$p.value, digits = 2),
                   '',
                   '',
                   '')
      }
      # metricnm <- rep(input$variablerecmet,5)
      group <- c('WTNE',
                 'TGNE',
                 'WTEE',
                 'TGEE',
                 '',
                 '',
                 ''
      )
      if(input$variablecomp=="Age effect"){
        group <- c('Wild-Type_2mo',
                   'Wild-Type_6mo',
                   'TgDyrk1A_2mo',
                   'TgDyrk1A_6mo',
                   '',
                   '',
                   ''
        )
      }
      means <- c(
        round(mean(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      sds <- c(
        round(sd(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      
      contrast <- c(
        'Genotype',
        'Treatment',
        'Genotype:Treatment',
        'WTNE - TGNE',
        'WTNE - WTEE',
        'WTNE - TGEE',
        'TGNE - TGEE'
      )
      if(input$variablecomp=="Age effect"){
        contrast <- c(
          'Genotype',
          'Age',
          'Genotype:Age',
          'Wild-Type_2mo - TgDyrk1A_2mo',
          'Wild-Type_2mo - Wild-Type_6mo',
          'Wild-Type_2mo - TgDyrk1A_6mo',
          'TgDyrk1A_2mo - TgDyrk1A_6mo'
        )
      }
      # anno <- cbind(metricnm,group,means,sds,pshap,anno)
      anno <- cbind(group,means,sds,pshap,contrast,annoprs)
      
      print(anno)
      # print(xtable(anno, caption = "Statistics for..."), include.rownames = FALSE)
      xtable(anno, caption = "Statistics for...")
    },spacing = 'm')
    
    #############################
    # Boxplots other metrics
    # output$metric <-  renderPlot({
    plotmetric <-  reactive({
      data = read.csv('Meri_metrics.csv')
      metname <- input$variablemet
      
      if(input$variablecomp2=="EE 2mo"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_2mo')
          data2 <- subset(data0,  data0$metric=='spine_C_2mo')
          data$value <- data$value+data2$value
        }
        else{
          data <- subset(data,data$metric==paste0(metname,'_2mo'))
        }
        data$gen <- NULL
        data$treat <- NULL
        data$gen[data$genotype=="WTNE"] <- 'WT'
        data$gen[data$genotype=="WTEE"] <- 'WT'
        data$gen[data$genotype=="TGNE"] <- 'TG'
        data$gen[data$genotype=="TGEE"] <- 'TG'
        data$treat[data$genotype=="WTNE"] <- 'NE'
        data$treat[data$genotype=="WTEE"] <- 'EE'
        data$treat[data$genotype=="TGNE"] <- 'NE'
        data$treat[data$genotype=="TGEE"] <- 'EE'
        
        data$genotype <- factor(data$genotype,levels=c('WTNE','TGNE','WTEE','TGEE'))
        levels(data$genotype)[levels(data$genotype)=="WTNE"] <- "WT"
        levels(data$genotype)[levels(data$genotype)=="TGNE"] <- "TgDyrk1A"
        levels(data$genotype)[levels(data$genotype)=="WTEE"] <- "WT_EE"
        levels(data$genotype)[levels(data$genotype)=="TGEE"] <- "TgDyrk1A_EE"
      }
      else if(input$variablecomp2=="EE 6mo"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_6mo')
          data2 <- subset(data0,  data0$metric=='spine_C_6mo')
          data$value <- data$value+data2$value
        }
        else{
          data <- subset(data,data$metric==paste0(metname,'_6mo'))
        }
        
        data$gen <- NULL
        data$treat <- NULL
        data$gen[data$genotype=="WTNE"] <- 'WT'
        data$gen[data$genotype=="WTEE"] <- 'WT'
        data$gen[data$genotype=="TGNE"] <- 'TG'
        data$gen[data$genotype=="TGEE"] <- 'TG'
        data$treat[data$genotype=="WTNE"] <- 'NE'
        data$treat[data$genotype=="WTEE"] <- 'EE'
        data$treat[data$genotype=="TGNE"] <- 'NE'
        data$treat[data$genotype=="TGEE"] <- 'EE'
        
        data$genotype <- factor(data$genotype,levels=c('WTNE','TGNE','WTEE','TGEE'))
        levels(data$genotype)[levels(data$genotype)=="WTNE"] <- "WT"
        levels(data$genotype)[levels(data$genotype)=="TGNE"] <- "TgDyrk1A"
        levels(data$genotype)[levels(data$genotype)=="WTEE"] <- "WT_EE"
        levels(data$genotype)[levels(data$genotype)=="TGEE"] <- "TgDyrk1A_EE"
      }
      else if(input$variablecomp2=="Age effect"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_2mo')
          data2 <- subset(data0,  data0$metric=='spine_C_2mo')
          data$value <- data$value+data2$value
          data1 <- data
          
          data <- subset(data0,  data0$metric=='spine_B_6mo')
          data2 <- subset(data0,  data0$metric=='spine_C_6mo')
          data$value <- data$value+data2$value
          data2 <- data
        }
        else{
          data1 <- subset(data,data$metric==paste0(metname,'_2mo'))
          data2 <- subset(data,data$metric==paste0(metname,'_6mo'))
        }
        
        data1$gen <- NULL
        data1$age <- NULL
        data1$treat <- NULL
        data1$gen[data1$genotype=="WTNE"] <- 'WT'
        data1$gen[data1$genotype=="WTEE"] <- 'WT'
        data1$gen[data1$genotype=="TGNE"] <- 'TG'
        data1$gen[data1$genotype=="TGEE"] <- 'TG'
        data1$treat[data1$genotype=="WTNE"] <- 'NE'
        data1$treat[data1$genotype=="WTEE"] <- 'EE'
        data1$treat[data1$genotype=="TGNE"] <- 'NE'
        data1$treat[data1$genotype=="TGEE"] <- 'EE'
        data1$age <- '2mo'
        
        data1$genotype <- as.factor(data1$genotype)
        # data1 <- subset(data1,data1$genotype=='WTNE' | data1$genotype=='TGNE') 
        levels(data1$genotype)[levels(data1$genotype)=="WTNE"] <- 'WTNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="TGNE"] <- 'TGNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="WTEE"] <- 'WTEE'
        levels(data1$genotype)[levels(data1$genotype)=="TGEE"] <- 'TGEE'
        
        data2$gen <- NULL
        data2$age <- NULL
        data2$treat <- NULL
        data2$gen[data2$genotype=="WTNE"] <- 'WT'
        data2$gen[data2$genotype=="WTEE"] <- 'WT'
        data2$gen[data2$genotype=="TGNE"] <- 'TG'
        data2$gen[data2$genotype=="TGEE"] <- 'TG'
        data2$treat[data2$genotype=="WTNE"] <- 'NE'
        data2$treat[data2$genotype=="WTEE"] <- 'EE'
        data2$treat[data2$genotype=="TGNE"] <- 'NE'
        data2$treat[data2$genotype=="TGEE"] <- 'EE'
        data2$age <- '6mo'
        data2$genotype <- as.factor(data2$genotype)
        # data2 <- subset(data2,data2$genotype=='WTNE' | data2$genotype=='TGNE') 
        levels(data2$genotype)[levels(data2$genotype)=="WTNE"] <- 'WTNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="TGNE"] <- 'TGNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="WTEE"] <- 'WTEEdis'
        levels(data2$genotype)[levels(data2$genotype)=="TGEE"] <- 'TGEEdis'
        
        data <- rbind(data1,data2)
        
        # data$genotype <- factor(data$genotype,levels=c('Wild-Type_2mo','Wild-Type_6mo','TgDyrk1A_2mo','TgDyrk1A_6mo'))
        data$genotype <- factor(data$genotype,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
        
        data$alphav <- NULL
        data$alphav[data$genotype %in% c('WTNE_2mo','TGNE_2mo','WTEE','TGEE')] <- 1
        data$alphav[data$genotype %in% c('WTNE_6mo','TGNE_6mo','WTEEdis','TGEEdis')] <- 0.6
      }
      
      data$metric <- data$value
      if(input$variablemet=='spine_A'|input$variablemet=='spine_B'|input$variablemet=='spine_C'|input$variablemet=='spine_B&C')
      {
        data$metric <- data$metric/20
        data$value <- data$value/20
      }
      data$blocking <- data$Animal
      data$group <- data$genotype
      df <- data
  
      df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      
      if(input$variablemet=='OR' | input$variablemet=='vglutvgat' | input$variablemet=='DYRK1A_act' | input$variablemet=='ratio_fEPSP'){
        if(input$variablecomp2!="Age effect"){
          df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking),gen=df$gen,treat=df$treat)
          print(summary(aov(values~gen*treat,data=df_st)))
          print(pairs(emmeans(aov(values~gen*treat,data=df_st),~gen*treat), adjust="holm"))
          anno <- as.data.frame(pairs(emmeans(aov(values~Genotypes,data=df_st), "Genotypes"), adjust = "holm"))$p.value[c(1,2,3,5)]
          }
        else{
          df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking),gen=df$gen,age=df$age,treat=df$treat)
          print(summary(aov(values~gen*age*treat,data=df_st)))
          print(as.data.frame(joint_tests(aov(values~gen*age*treat,data=df_st), by = "age")))
          print(as.data.frame(joint_tests(aov(values~gen*age*treat,data=df_st), by = "gen")))
          print(as.data.frame(joint_tests(aov(values~gen*age*treat,data=df_st), by = "treat")))
          print(as.data.frame(pairs(emmeans(aov(values~gen*age*treat,data=df_st), ~gen*age*treat), adjust = "holm")))
          anno <- as.data.frame(pairs(emmeans(aov(values~gen*age*treat,data=df_st), ~gen*age*treat), adjust = "holm"))$p.value
          
          # Age
          # print("Age")
          # print(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]))
          # print(summary(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),])))
          # print(pairs(emmeans(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]), ~gen*age), adjust = "holm"))
          # # Discontinuation
          # print("Discontinuation")
          # print(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]))
          # print(summary(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),])))
          # print(pairs(emmeans(aov(values~gen*age,data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]), ~gen*age), adjust = "holm"))
          
        }
      }
      else{
        if(input$variablecomp2!="Age effect"){
          summary(mixed(values~Genotypes+(1|blocking),data=df_st))
          anno <- as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), ~Genotypes), adjust = "holm"))$p.value[c(1,2,3,5)]
        }
        else{
          df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking),gen=df$gen,age=df$age,treat=df$treat)
          print(summary(mixed(values~gen*age*treat+(1|blocking),data=df_st)))
          print(as.data.frame(joint_tests(mixed(values~gen*age*treat+(1|blocking),data=df_st), by = "age")))
          print(as.data.frame(joint_tests(mixed(values~gen*age*treat+(1|blocking),data=df_st), by = "gen")))
          print(as.data.frame(joint_tests(mixed(values~gen*age*treat+(1|blocking),data=df_st), by = "treat")))
          
          print(as.data.frame(pairs(emmeans(mixed(values~gen*age*treat+(1|blocking),data=df_st), ~gen*age*treat), adjust = "holm")))
          anno <- as.data.frame(pairs(emmeans(mixed(values~gen*age*treat+(1|blocking),data=df_st), ~gen*age*treat), adjust = "holm"))$p.value
          
          # Age
          # print("Age")
          # print(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]))
          # print(summary(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),])))
          # print(pairs(emmeans(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]), ~gen*age), adjust = "holm"))
          # # Discontinuation
          # print("Discontinuation")
          # print(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]))
          # print(summary(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),])))
          # print(pairs(emmeans(mixed(values~gen*age+(1|blocking),data=df_st[df_st$Genotypes %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]), ~gen*age), adjust = "holm"))
        }
      }
        # ggviolin(df, x="group", y="metric", fill = "group",width=0.5, 
        #          trim=TRUE,
        #          palette="Set1",
        #          add = c("boxplot","dotplot"),
        #          add.params = list(size=0.5),
        #          ylab = metname) +
        # anno <- as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), "Genotypes"), adjust = "fdr"))$p.value[c(1,2,3,5)]
        if(input$variablecomp2!="Age effect"){
          ggbarplot(df, x="group", y="metric", fill = "group", color = "group",width=0.3, 
                    trim=TRUE,
                    palette="Set1",
                    # add = c("boxplot","dotplot"),
                    add = c("mean_se"),
                    add.params = list(size=0.5),
                    ylab = metname) +
            geom_hline(yintercept=mean(df$metric[df$group=="WT"]),linetype=2) +
            # stat_signif(annotation='',#formatC(paste0(groupanno3[1],groupanno3[2],groupanno3[3])),
            #             textsize=6,
            #             tip_length = 0,
            #             # position = position_identity(),
            #             margin_top = 0.7,
            #             step_increase = 0.15,
            #             # y_position = c((mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.9,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.8,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.7),
            #             y_position = (mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.65,
            #             comparisons=list(c(1,4))
            # ) +
            theme(axis.text=element_text(size=20),
                  axis.title = element_text(size=25),
                  title=element_text(size=25),
                  legend.text = element_text(size=25)) +
            scale_y_continuous(expand = expand_scale(mult = c(0, .1)),breaks = scales::pretty_breaks(n = 6)) #+
            # geom_signif(annotation=formatC(anno, digits=2),
            #             aes(group = group), 
            #                    # test="t.test",
            #                    # test=as.list(as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df), "Genotypes"), adjust = "fdr"))),
            #                    # test=pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), "Genotypes"), adjust = "fdr"),
            #                    # label = "p.signif", 
            #                    map_signif_level = T,
            #                    textsize=6,
            #                   y_position = c((mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.1,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.25,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.4,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.55),
            #                   comparisons=list(c('WT','TgDyrk1A'),
            #                                    c('WT','WT_EE'),
            #                                    c('WT','TgDyrk1A_EE'),
            #                                    c('TgDyrk1A','TgDyrk1A_EE')
            #                   )
            #                   # comparisons=list(c('WTNE','TGNE'),
            #                   #                  c('WTNE','WTEE'),
            #                   #                  c('WTNE','TGEE'),
            #                   #                  c('TGNE','TGEE')
            #                   # )
            # )
        }
        else{
          palettev <- c("#E41A1C","#E41A1C","#377EB8","#377EB8","#4DAF4A","#4DAF4A","#984EA3","#984EA3")
          print(df$alphav)
          ggbarplot(df, x="group", y="metric", fill = "group", color = "group",width=0.3, 
                    trim=TRUE,
                    # palette="8-class Paired",
                    palette=palettev,
                    # alpha="alphav",
                    # add = c("boxplot","dotplot"),
                    add = c("mean_se"),
                    add.params = list(size=0.5),
                    ylab = metname) +
            geom_hline(yintercept=mean(df$metric[df$group=="WTNE_2mo"]),linetype=2) +
            # stat_signif(annotation='',#formatC(paste0(groupanno3[1],groupanno3[2],groupanno3[3])),
            #             textsize=6,
            #             tip_length = 0,
            #             # position = position_identity(),
            #             margin_top = 0.7,
            #             step_increase = 0.15,
            #             # y_position = c((mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.9,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.8,(mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.7),
            #             y_position = (mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.65,
            #             comparisons=list(c(1,8))
            # ) +
            theme(axis.text=element_text(size=20),
                  axis.title = element_text(size=25),
                  title=element_text(size=25),
                  legend.text = element_text(size=25)) +
            scale_y_continuous(expand = expand_scale(mult = c(0, .1)),breaks = scales::pretty_breaks(n = 6)) #+
            # geom_signif(annotation=formatC(anno[c(27,24,9,2)], digits=2),
            #             aes(group = group), 
            #             # test="t.test",
            #             # test=as.list(as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df), "Genotypes"), adjust = "fdr"))),
            #             # test=pairs(emmeans(mixed(values~gen*age*treat+(1|blocking),data=df_st), ~gen*age*treat), adjust = "holm"),
            #             # label = "p.signif", 
            #             map_signif_level = T,
            #             textsize=6,
            #             y_position = rep((mean(df$metric[df$genotype==df$genotype[df$metric==max(df$metric)]])+std.error(df$metric))*1.4,4),
            #             comparisons=list(c(1,2),
            #                              c(3,4),
            #                              c(5,6),
            #                              c(7,8)
            #             )
            #             # comparisons=list(c('WTNE','TGNE'),
            #             #                  c('WTNE','WTEE'),
            #             #                  c('WTNE','TGEE'),
            #             #                  c('TGNE','TGEE')
            #             # )
            # )
        }
      
      # recordedplot <- recordPlot()
      # save(recordedplot,file='./templot.Rdata')
      
      # }
      # else if(input$variablesubl == 'Stratum Radiatum'){
      #   df$metric <- df$str_rad
      #   bvdplot(df,'Stratum Radiatum width')
      # }
      # else if(input$variablesubl == 'Stratum Lacunosum'){
      #   df$metric <- df$str_lac
      #   bvdplot(df,'Stratum Lacunosum width')
      # }
      #   )
      # )
    })
    
    output$metric <-  renderPlot({
      if(input$variablemet=='spine_A'|input$variablemet=='spine_B'|input$variablemet=='spine_C')
      {
        print(plotmetric()+scale_y_continuous(expand = expand_scale(mult = c(0, .1)),limits=c(0,1.45),breaks=seq(0,1.2,by=0.2)))
      # }
      # else if(input$variablemet=='ratio_fEPSP'){
      #   print(plotmetric()+scale_y_continuous(expand = expand_scale(add = 1),limits=c(1,2),breaks=seq(1,2,by=0.2)))
      }else{
        print(plotmetric())
      }
    })
    
    output$downloadPlot <- downloadHandler(
      filename = function(){paste(input$variablemet,'_',input$variablecomp2, '.pdf', sep = '')},
      
      content = function(file){
        if(input$variablemet=='spine_A'|input$variablemet=='spine_B'|input$variablemet=='spine_C')
        {
          ggsave(file,width=3,height=8,plotmetric()+scale_y_continuous(expand = expand_scale(mult = c(0, .1)),limits=c(0,1.45),breaks=seq(0,1.2,by=0.2)),useDingbats=F)
        }
        else if(input$variablecomp2=="Age effect")
        {
          ggsave(file,width=12,height=6,plotmetric(),useDingbats=F)
        }
        else{
          ggsave(file,width=6,height=8,plotmetric(),useDingbats=F)
        }
        },
      
      contentType = "application/pdf"
    )
      
    output$table <- renderTable({
      data = read.csv('Meri_metrics.csv')
      metname <- input$variablemet
      
      if(input$variablecomp2=="EE 2mo"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_2mo')
          data2 <- subset(data0,  data0$metric=='spine_C_2mo')
          data$value <- data$value+data2$value
          data1 <- data
        }
        else{
          data <- subset(data,data$metric==paste0(metname,'_2mo'))
        }
        data$genotype <- factor(data$genotype,levels=c('WTNE','TGNE','WTEE','TGEE'))
        # levels(data$genotype)[levels(data$genotype)=="WTNE"] <- "WT"
        # levels(data$genotype)[levels(data$genotype)=="TGNE"] <- "TgDyrk1A"
        # levels(data$genotype)[levels(data$genotype)=="WTEE"] <- "WT_EE"
        # levels(data$genotype)[levels(data$genotype)=="TGEE"] <- "TgDyrk1A_EE"
      }
      else if(input$variablecomp2=="EE 6mo"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_6mo')
          data2 <- subset(data0,  data0$metric=='spine_C_6mo')
          data$value <- data$value+data2$value
          data1 <- data
        }
        else{
          data <- subset(data,data$metric==paste0(metname,'_6mo'))
        }
        data$genotype <- factor(data$genotype,levels=c('WTNE','TGNE','WTEE','TGEE'))
        # levels(data$genotype)[levels(data$genotype)=="WTNE"] <- "WT"
        # levels(data$genotype)[levels(data$genotype)=="TGNE"] <- "TgDyrk1A"
        # levels(data$genotype)[levels(data$genotype)=="WTEE"] <- "WT_EE"
        # levels(data$genotype)[levels(data$genotype)=="TGEE"] <- "TgDyrk1A_EE"
      }
      else if(input$variablecomp2=="Age effect"){
        if(metname == 'spine_B&C'){
          data0 <- data
          data <- subset(data,  data$metric=='spine_B_2mo')
          data2 <- subset(data0,  data0$metric=='spine_C_2mo')
          data$value <- data$value+data2$value
          data1 <- data
          
          data <- subset(data0,  data0$metric=='spine_B_6mo')
          data2 <- subset(data0,  data0$metric=='spine_C_6mo')
          data$value <- data$value+data2$value
          data2 <- data
        }
        else{
          data1 <- subset(data,data$metric==paste0(metname,'_2mo'))
          data2 <- subset(data,data$metric==paste0(metname,'_6mo'))
        }
        
        data1 <- subset(data1,data1$genotype=='WTNE' | data1$genotype=='TGNE') 
        levels(data1$genotype)[levels(data1$genotype)=="WTNE"] <- 'Wild-Type_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="TGNE"] <- 'TgDyrk1A_2mo'
        data2 <- subset(data2,data2$genotype=='WTNE' | data2$genotype=='TGNE') 
        levels(data2$genotype)[levels(data2$genotype)=="WTNE"] <- 'Wild-Type_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="TGNE"] <- 'TgDyrk1A_6mo'
        
        data <- rbind(data1,data2)
        
        data$genotype <- factor(data$genotype,levels=c('Wild-Type_2mo','Wild-Type_6mo','TgDyrk1A_2mo','TgDyrk1A_6mo'))
        # lapply(strsplit(string, ":"), `[`, 2)
        data$gen <- sapply(strsplit(as.character(data$genotype),'_'),'[',1)
        data$age <- sapply(strsplit(as.character(data$genotype),'_'),'[',2)
      }
      
      data$metric <- data$value
      if(input$variablemet=='spine_A'|input$variablemet=='spine_B'|input$variablemet=='spine_C'|input$variablemet=='spine_B&C')
      {
        data$metric <- data$metric/20
        data$value <- data$value/20
      }
      data$blocking <- data$Animal
      data$group <- data$genotype
      df <- data
      
      if(input$variablecomp2=="Age effect"){
        df$group <- factor(df$group,levels=c('Wild-Type_2mo','Wild-Type_6mo','TgDyrk1A_2mo','TgDyrk1A_6mo'))
      }
      else{
        # df$group <- factor(df$group,levels=c('WT','TgDyrk1A','WT_EE','TgDyrk1A_EE'))
        df$gen <- sapply(df$group,substring,first=1,last=2)
        df$treat <- sapply(df$group,substring,first=3,last=4)
      }
      
      df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      if(input$variablemet=='OR' | input$variablemet=='vglutvgat' | input$variablemet=='DYRK1A_act' | input$variablemet=='ratio_fEPSP'){
        if(input$variablecomp2=="Age effect"){
          print(summary(aov(metric~gen*age,data=df)))
          groupannodf <- as.data.frame(summary(aov(metric~gen*age,data=df))[[1]])[1:3,1]
          groupannoF <- as.data.frame(summary(aov(metric~gen*age,data=df))[[1]])[1:3,4]
          groupannop <- as.data.frame(summary(aov(metric~gen*age,data=df))[[1]])[1:3,5]
        }
        else{
          print(summary(aov(metric~gen*treat,data=df)))
          groupannodf <- as.data.frame(summary(aov(metric~gen*treat,data=df))[[1]])[1:3,1]
          groupannoF <- as.data.frame(summary(aov(metric~gen*treat,data=df))[[1]])[1:3,4]
          groupannop <- as.data.frame(summary(aov(metric~gen*treat,data=df))[[1]])[1:3,5]
        }
      }
      else{
        if(input$variablecomp2=="Age effect"){
          print(summary(mixed(metric~gen*age+(1|blocking),data=df)))
          groupannodf <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.den.Df
          groupannoF <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.F
          groupannop <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.Pr..F.
        }
        else{
          # print(summary(mixed(metric~gen*treat+(1|blocking),data=df)))
          groupannodf <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.den.Df
          groupannoF <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.F
          groupannop <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
          print(as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1]))
        }
      }
      
      # Kruskal test to take then format
      stat.all <- compare_means(
        metric ~ group, data = df,
        method = "kruskal.test"
      )
      
      # Pairwise wilcox-test between groups to take then format
      stat.pairwise <- compare_means(
        metric ~ group, data = df,
        method = "wilcox.test"
      ) %>%
        mutate(y.position = max(df$metric)*seq(1.05,1.55,length.out=6))
      
      stat.pairwise$y.position[6] <- stat.pairwise$y.position[4]
      stat.pairwise$y.position[5] <- stat.pairwise$y.position[4]
      stat.pairwise <- stat.pairwise[c(1:3,5),]
      if(input$variablecomp=="Age effect"){
        stat.pairwise <- stat.pairwise[c(1:3,6),]
      }
      
      # Linear mixed effects model
      df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      # m1 <- aov(value~gen*treat,data=df)
      # print(emm1 <- emmeans(m1, ~gen*treat))
      if(input$variablemet=='OR' | input$variablemet=='vglutvgat' | input$variablemet=='DYRK1A_act' | input$variablemet=='ratio_fEPSP'){
        m1 <- aov(values~Genotypes,data=df_st)
        (emm1 <- emmeans(m1, "Genotypes"))
        prs <- as.data.frame(pairs(emm1, adjust = "holm"))
        print(prs)
      }
      else{
        m1 <- mixed(values~Genotypes+(1|blocking),data=df_st)
        (emm1 <- emmeans(m1, ~Genotypes))
        prs <- as.data.frame(pairs(emm1, adjust = "holm"))
        print(prs)
      }
        
      #coefs
      shap_wtne <- shapiro.test(df$metric[df$group==unique(df$group)[1]])
      shap_tgne <- shapiro.test(df$metric[df$group==unique(df$group)[2]])
      shap_wtee <- shapiro.test(df$metric[df$group==unique(df$group)[3]])
      shap_tgee <- shapiro.test(df$metric[df$group==unique(df$group)[4]])
      annoprs0 <- prs[c(1:3,5),c(1,4,5,6)]
      if(input$variablecomp2=="Age effect"){
        annoprs0 <- prs[c(1:3,6),c(1,4,5,6)]
      }
      annoprs <- data.frame()[1:7,]
      annoprs$df <- c(formatC(groupannodf,digits=2),round(annoprs0$df, digits=2))
      annoprs$t.ratio <- c(formatC(groupannoF,digits=2),round(annoprs0$t.ratio, digits=2))
      annoprs$p.value <- c(formatC(groupannop,digits=2),formatC(annoprs0$p.value,digits=2,format = 'g'))
      
      annoprs <- rename(as.data.frame(annoprs),'F/t.ratio'='t.ratio')
      pshap <- c(round(shap_wtne$p.value, digits = 2),
                 round(shap_tgne$p.value, digits = 2),
                 round(shap_wtee$p.value, digits = 2),
                 round(shap_tgee$p.value, digits = 2),
                 '',
                 '',
                 '')
      # metricnm <- rep(input$variablerecmet,5)
      group <- c('WTNE',
                 'TGNE',
                 'WTEE',
                 'TGEE',
                 '',
                 '',
                 ''
      )
      if(input$variablecomp2=="Age effect"){
        group <- c('Wild-Type_2mo',
                   'Wild-Type_6mo',
                   'TgDyrk1A_2mo',
                   'TgDyrk1A_6mo',
                   '',
                   '',
                   ''
        )
      }
      means <- c(
        round(mean(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      sds <- c(
        round(sd(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      if(input$variablecomp2=="Age effect"){
        contrast <- c(
          'Genotype',
          'Age',
          'Genotype:Age',
          'Wild-Type_2mo - Wild-Type_6mo',
          'Wild-Type_2mo - TgDyrk1A_2mo',
          'Wild-Type_2mo - TgDyrk1A_6mo',
          'TgDyrk1A_2mo - TgDyrk1A_6mo'
        )
      }
      else{
        contrast <- c(
          'Genotype',
          'Treatment',
          'Genotype:Treatment',
          'WTNE - TGNE',
          'WTNE - WTEE',
          'WTNE - TGEE',
          'TGNE - TGEE'
        )
      }
      # anno <- cbind(metricnm,group,means,sds,pshap,anno)
      anno <- cbind(group,means,sds,pshap,contrast,annoprs)
      
      print(anno)
      # print(xtable(anno, caption = "Statistics for..."), include.rownames = FALSE)
      xtable(anno, caption = "Statistics for...")
    },spacing = 'm')
    
    
    
    
    # output$recmet <- renderPlot({
    plotrecmetric <-  reactive({
      if(input$variablecomp3 == 'EE 2mo'){
        data <- read.csv('nm_data_ap_ee_0.7bound.csv')
        data$genotype <- as.factor(data$genotype)
        data$genotype <- factor(data$genotype,levels=c(levels(data$genotype),c('WTNE','TGNE','WTEE','TGEE')))
        levels(data$genotype)[levels(data$genotype)=="Wild-Type"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_EE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_EE"] <- "TGEE"
      }
      else if(input$variablecomp3 == 'EE 6mo'){
        data <- read.csv('nm_data_ap_ee_0.7bound_withdisc.csv')
        data$genotype <- as.factor(data$genotype)
        data$genotype <- factor(data$genotype,levels=c(levels(data$genotype),c('WTNE','TGNE','WTEE','TGEE')))
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_NE_NE"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_NE_NE"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_EE_NE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_EE_NE"] <- "TGEE"
      }
      else if(input$variablecomp3 == 'Age effect'){
        data1 <- read.csv('nm_data_ap_ee_0.7bound.csv')
        data2 <- read.csv('nm_data_ap_ee_0.7bound_withdisc.csv')
        
        data1$gen <- NULL
        data1$age <- NULL
        data1$treat <- NULL
        data1$gen[data1$genotype=="Wild-Type"] <- 'WT'
        data1$gen[data1$genotype=="Wild-Type_EE"] <- 'WT'
        data1$gen[data1$genotype=="TgDyrk1A"] <- 'TG'
        data1$gen[data1$genotype=="TgDyrk1A_EE"] <- 'TG'
        data1$treat[data1$genotype=="Wild-Type"] <- 'NE'
        data1$treat[data1$genotype=="Wild-Type_EE"] <- 'EE'
        data1$treat[data1$genotype=="TgDyrk1A"] <- 'NE'
        data1$treat[data1$genotype=="TgDyrk1A_EE"] <- 'EE'
        data1$age <- '2mo'
        
        data1$genotype <- as.factor(data1$genotype)
        # data1 <- subset(data1,data1$genotype=='WTNE' | data1$genotype=='TGNE') 
        levels(data1$genotype)[levels(data1$genotype)=="Wild-Type"] <- 'WTNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="TgDyrk1A"] <- 'TGNE_2mo'
        levels(data1$genotype)[levels(data1$genotype)=="Wild-Type_EE"] <- 'WTEE'
        levels(data1$genotype)[levels(data1$genotype)=="TgDyrk1A_EE"] <- 'TGEE'
        
        data2$gen <- NULL
        data2$age <- NULL
        data2$treat <- NULL
        data2$gen[data2$genotype=="Wild-Type_NE_NE"] <- 'WT'
        data2$gen[data2$genotype=="Wild-Type_EE_NE"] <- 'WT'
        data2$gen[data2$genotype=="TgDyrk1A_NE_NE"] <- 'TG'
        data2$gen[data2$genotype=="TgDyrk1A_EE_NE"] <- 'TG'
        data2$treat[data2$genotype=="Wild-Type_NE_NE"] <- 'NE'
        data2$treat[data2$genotype=="Wild-Type_EE_NE"] <- 'EE'
        data2$treat[data2$genotype=="TgDyrk1A_NE_NE"] <- 'NE'
        data2$treat[data2$genotype=="TgDyrk1A_EE_NE"] <- 'EE'
        data2$age <- '6mo'
        data2$genotype <- as.factor(data2$genotype)
        # data2 <- subset(data2,data2$genotype=='WTNE' | data2$genotype=='TGNE') 
        levels(data2$genotype)[levels(data2$genotype)=="Wild-Type_NE_NE"] <- 'WTNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="TgDyrk1A_NE_NE"] <- 'TGNE_6mo'
        levels(data2$genotype)[levels(data2$genotype)=="Wild-Type_EE_NE"] <- 'WTEEdis'
        levels(data2$genotype)[levels(data2$genotype)=="TgDyrk1A_EE_NE"] <- 'TGEEdis'
        
        data=rbind(data1,data2)
        data$genotype <- factor(data$genotype,levels=c('WTNE_2mo','WTNE_6mo','TGNE_2mo','TGNE_6mo','WTEE','WTEEdis','TGEE','TGEEdis'))
        
        data$alphav <- NULL
        data$alphav[data$genotype %in% c('WTNE_2mo','TGNE_2mo','WTEE','TGEE')] <- 1
        data$alphav[data$genotype %in% c('WTNE_6mo','TGNE_6mo','WTEEdis','TGEEdis')] <- 0.6
      }
      
      # df <- data
      # df$gen <- sapply(df$genotype,substring,first=1,last=2)
      # df$treat <- sapply(df$genotype,substring,first=3,last=4)
      # df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,treat=df$treat)
      # # df$group <- factor(df$group,levels=c('WT','TgDyrk1A','WT_EE','TgDyrk1A_EE'))
      if(input$variablecomp3 != 'Age effect'){
        df <- data
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$treat <- sapply(df$genotype,substring,first=3,last=4)
        df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,treat=df$treat)
        df$group <- factor(df$group,levels=c('WTNE','TGNE','WTEE','TGEE'))
        df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*treat+(1|blocking),data=df))
        print(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))
        groupanno <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Treatment effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Treatment interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(6,5,3,2)]
        }
      else{
        df <- data
        # df$gen <- sapply(df$genotype,substring,first=1,last=2)
        # df$age <- sapply(df$genotype,substring,first=4,last=6)
        df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,age=df$age,treat=df$treat,alphav=df$alphav)
        
        # df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*age*treat+(1|blocking),data=df))
        print(summary(mixed(metric~gen*age*treat+(1|blocking),data=df)))
        print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "age")))
        print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "gen")))
        print(as.data.frame(joint_tests(mixed(metric~gen*age*treat+(1|blocking),data=df), by = "treat")))
        print(as.data.frame(pairs(emmeans(mixed(metric~gen*age*treat+(1|blocking),data=df), ~gen*age*treat), adjust = "holm")))
        
        # Age
        #print("Age")
        #print(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]))
        # print(summary(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),])))
        # print(pairs(emmeans(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTNE_2mo","WTNE_6mo","TGNE_2mo","TGNE_6mo"),]), ~gen*age), adjust = "holm"))
        # Discontinuation
        # print("Discontinuation")
        # print(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]))
        # print(summary(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),])))
        # print(pairs(emmeans(mixed(metric~gen*age+(1|blocking),data=df[df$group %in% c("WTEE","WTEEdis","TGEE","TGEEdis"),]), ~gen*age), adjust = "holm"))
        
        groupanno <- as.data.frame(mixed(metric~gen*age*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Age effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Age interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*age*treat+(1|blocking),data=df), ~gen*age*treat), adjust = "holm"))$p.value
      }
      
      # groupanno3 <- groupanno
      # groupanno3[1] <- as.character(symnum(as.numeric(groupanno[1]), corr = FALSE,
      #                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                     symbols = c("????????", "??????", "????", "??", "")))
      # groupanno3[2] <- as.character(symnum(as.numeric(groupanno[2]), corr = FALSE,
      #                                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                      symbols = c("????????", "??????", "????", "??", "")))
      # groupanno3[3] <- as.character(symnum(as.numeric(groupanno[3]), corr = FALSE,
      #                                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      #                                      symbols = c("????????", "??????", "????", "??", "")))
      # print(paste0(groupanno3[1],groupanno3[2],groupanno3[3]))
      # anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(1,2,3,5)]
      # anno <- as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), "Genotypes"), adjust = "holm"))$p.value[c(1,2,3,5)]
      # anno <- as.data.frame(pairs(emmeans(mixed(values~Genotypes,data=df_st), "Genotypes"), adjust = "holm"))$p.value[c(1,2,3,5)]
      
      
      if(input$variablecomp3 != 'Age effect'){
        pl <- ggviolin(df, x="group", y="metric", fill = "group",width=0.3, 
                       trim=TRUE,
                       palette="Set1",
                       add = c("boxplot","dotplot"),
                       # add = c("mean_se"),
                       add.params = list(size=0.5),
                       ylab = input$variablerecmet) +
          geom_hline(yintercept=median(df$metric[df$group=="WTNE"]),linetype=2) +
          theme(axis.text=element_text(size=20),
                axis.title = element_text(size=25),
                title=element_text(size=25),
                legend.text = element_text(size=25)) #+
          # stat_signif(annotation='',#formatC(paste0(groupanno3[1],groupanno3[2],groupanno3[3])),
          #             textsize=6,
          #             tip_length = 0,
          #             # position = position_identity(),
          #             margin_top = 0.7,
          #             step_increase = 0.15,
          #             # y_position = c(max(df$metric)*1.9,max(df$metric)*1.8,max(df$metric)*1.7),
          #             y_position = max(df$metric)*1.5,
          #             comparisons=list(c(1,4))
          # ) +
          
          # annotate("text", label = groupanno2[4], x=1, y=max(df$metric)*1.55, size=6, hjust=0) +
          # annotate("text", label = groupanno2[3], x=1, y=max(df$metric)*1.65, size=6, hjust=0) +
          # annotate("text", label = groupanno2[2], x=1, y=max(df$metric)*1.75, size=6, hjust=0) +
          # annotate("text", label = groupanno2[1], x=1, y=max(df$metric)*1.85, size=6, hjust=0)
        
        pl
        # pl + geom_signif(annotation=formatC(anno, digits=2),
        #                aes(group = group), 
        #                # test="t.test",
        #                # test=as.list(as.data.frame(pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df), "Genotypes"), adjust = "fdr"))),
        #                # test=pairs(emmeans(mixed(values~Genotypes+(1|blocking),data=df_st), "Genotypes"), adjust = "holm"),
        #                # label = "p.signif", 
        #                map_signif_level = T,
        #                textsize=6,
        #                # margin_top = 0.12,
        #                # step_increase = 0.1,
        #                y_position = c(max(df$metric)*1.05,max(df$metric)*1.15,max(df$metric)*1.25,max(df$metric)*1.35),
        #                # comparisons=list(c('WT','TgDyrk1A'),
        #                #                  c('WT','WT_EE'),
        #                #                  c('WT','TgDyrk1A_EE'),
        #                #                  c('TgDyrk1A','TgDyrk1A_EE')
        #                # )
        #                comparisons=list(c(1,2),
        #                                 c(1,3),
        #                                 c(1,4),
        #                                 c(2,4)
        #                ) 
        # )
      } else {
        palettev <- c("#E41A1C","#E41A1C","#377EB8","#377EB8","#4DAF4A","#4DAF4A","#984EA3","#984EA3")
        pl <- ggviolin(df, x="group", y="metric", fill = "group",width=0.3, 
                       trim=TRUE,
                       # palette="8-class Paired",
                       palette=palettev,
                       alpha="alphav",
                       add = c("boxplot","dotplot"),
                       # add = c("mean_se"),
                       add.params = list(size=0.5,alpha = "alphav"),
                       ylab = input$variablerecmet) #+
          # geom_hline(yintercept=median(df$metric[df$group=="WTNE"]),linetype=2) +
          # theme(axis.text=element_text(size=20),
          #       axis.title = element_text(size=25),
          #       title=element_text(size=25),
          #       legend.text = element_text(size=25)) #+
          # stat_signif(annotation='',#formatC(paste0(groupanno3[1],groupanno3[2],groupanno3[3])),
          #             textsize=6,
          #             tip_length = 0,
          #             # position = position_identity(),
          #             margin_top = 0.7,
          #             step_increase = 0.15,
          #             # y_position = c(max(df$metric)*1.9,max(df$metric)*1.8,max(df$metric)*1.7),
          #             y_position = max(df$metric)*1.5,
          #             comparisons=list(c(1,4))
          # ) +
          
          # annotate("text", label = groupanno2[4], x=1, y=max(df$metric)*1.55, size=6, hjust=0) +
          # annotate("text", label = groupanno2[3], x=1, y=max(df$metric)*1.65, size=6, hjust=0) +
          # annotate("text", label = groupanno2[2], x=1, y=max(df$metric)*1.75, size=6, hjust=0) +
          # annotate("text", label = groupanno2[1], x=1, y=max(df$metric)*1.85, size=6, hjust=0)
        
        pl
        # pl + geom_signif(annotation=formatC(anno[c(27,24,9,2)], digits=2),
        #                  aes(group = group), 
        #                  map_signif_level = T,
        #                  textsize=6,
        #                  # margin_top = 0.12,
        #                  # step_increase = 0.1,
        #                  y_position = rep(max(df$metric)*1.1,4),
        #                  comparisons=list(c(1,2),
        #                                   c(3,4),
        #                                   c(5,6),
        #                                   c(7,8)
        #                  ) 
        #                  
        # ) 
          # scale_fill_brewer(palette="8-class Paired")
      }
      # p$layers <- c(#stat_function(fun = function(x) x, mapping = aes(xmin=0,xmax=..x..), color='gray',linetype="dashed"),
      #   annotate("text", label = groupanno[3], x=1, y=max(df$metric)*1.5, size=6, hjust=0),
      #   annotate("text", label = groupanno[2], x=1, y=max(df$metric)*1.6, size=6, hjust=0),
      #   annotate("text", label = groupanno[1], x=1, y=max(df$metric)*1.7, size=6, hjust=0),
      #   p$layers)
      # p
    })
    
    output$recmet <-  renderPlot({
      print(plotrecmetric())
    })
    
    output$downloadPlot2 <- downloadHandler(
      filename = function(){paste(input$variablerecmet,'_',input$variablecomp3, '.pdf', sep = '')},
      
      content = function(file){
        ggsave(file,width=12,height=6,plotrecmetric(),useDingbats=F)
      },
      
      contentType = "application/pdf"
    )
    
    output$tablerec <- renderTable({
      if(input$variablecomp3 == 'EE 2mo'){
        data <- read.csv('nm_data_ap_ee_0.7bound.csv')
        levels(data$genotype)[levels(data$genotype)=="Wild-Type"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_EE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_EE"] <- "TGEE"
      }
      else if(input$variablecomp3 == 'EE 6mo'){
        data <- read.csv('nm_data_ap_ee_0.7bound_withdisc.csv')
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_NE_NE"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_NE_NE"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_EE_NE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_EE_NE"] <- "TGEE"
      }
      else if(input$variablecomp3 == 'Age effect'){
        data1 <- read.csv('nm_data_ap_ee_0.7bound.csv')
        data2 <- read.csv('nm_data_ap_ee_0.7bound_withdisc.csv')
        data=rbind(data1,data2)
        data <- subset(data,data$genotype!="Wild-Type_EE"&data$genotype!="TgDyrk1A_EE"&data$genotype!="Wild-Type_EE_NE"&data$genotype!="TgDyrk1A_EE_NE")
        levels(data$genotype)[levels(data$genotype)=="Wild-Type"] <- "WTNE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A"] <- "TGNE"
        levels(data$genotype)[levels(data$genotype)=="Wild-Type_NE_NE"] <- "WTEE"
        levels(data$genotype)[levels(data$genotype)=="TgDyrk1A_NE_NE"] <- "TGEE"
      }
      
      # df <- data
      # df$gen <- sapply(df$genotype,substring,first=1,last=2)
      # df$treat <- sapply(df$genotype,substring,first=3,last=4)
      # df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,treat=df$treat)
      # # df$group <- factor(df$group,levels=c('WT','TgDyrk1A','WT_EE','TgDyrk1A_EE'))
      # df$group <- factor(df$group,levels=c('WTNE','TGNE','WTEE','TGEE'))
      # df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      # 
      if(input$variablecomp3 != 'Age effect'){
        df <- data
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$treat <- sapply(df$genotype,substring,first=3,last=4)
        df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,treat=df$treat)
        df$group <- factor(df$group,levels=c('WTNE','TGNE','WTEE','TGEE'))
        df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*treat+(1|blocking),data=df))
        groupanno <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Treatment effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Treatment interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(6,5,3,2)]
        groupannodf <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.den.Df
        groupannoF <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.F
        groupannop <- as.data.frame(mixed(metric~gen*treat+(1|blocking),data=df)[1])$anova_table.Pr..F.
      }
      else{
        df <- data
        df$gen <- sapply(df$genotype,substring,first=1,last=2)
        df$age <- sapply(df$genotype,substring,first=3,last=4)
        df <- data.frame(metric=data[,which(names(data)==input$variablerecmet)],group=data$genotype,blocking=data$animal,gen=df$gen,age=df$age)
        df$group <- factor(df$group,levels=c('WTNE','WTEE','TGNE','TGEE'))
        levels(df$group)[levels(df$group)=="WTNE"] <- "Wild-Type_2mo"
        levels(df$group)[levels(df$group)=="TGNE"] <- "TgDyrk1A_2mo"
        levels(df$group)[levels(df$group)=="WTEE"] <- "Wild-Type_6mo"
        levels(df$group)[levels(df$group)=="TGEE"] <- "TgDyrk1A_6mo"
        df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
        print(mixed(metric~gen*age+(1|blocking),data=df))
        groupanno <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.Pr..F.
        groupanno2 <- groupanno
        groupanno2[1]<-'Mixed-effects KR approx.'
        groupanno2[2]<-paste0('Genotype effect p = ',formatC(groupanno[1], digits=2))
        groupanno2[3]<-paste0('Age effect p = ',formatC(as.numeric(groupanno[2]), digits=2))
        groupanno2[4]<-paste0('Genotype:Age interaction p = ',formatC(as.numeric(groupanno[3]), digits=2))
        # print(as.data.frame(pairs(emmeans(mixed(metric~gen*age+(1|blocking),data=df), ~gen*age), adjust = "holm"))$p.value)
        anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*age+(1|blocking),data=df), ~gen*age), adjust = "holm"))$p.value[c(5,2)]
        groupannodf <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.den.Df
        groupannoF <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.F
        groupannop <- as.data.frame(mixed(metric~gen*age+(1|blocking),data=df)[1])$anova_table.Pr..F.
      }
      # anno <- as.data.frame(pairs(emmeans(mixed(metric~gen*treat+(1|blocking),data=df), ~gen*treat), adjust = "holm"))$p.value[c(6,5,3,2)]
      
      
      
      # Kruskal test to take then format
      stat.all <- compare_means(
        metric ~ group, data = df,
        method = "kruskal.test"
      )
      
      # Pairwise wilcox-test between groups to take then format
      stat.pairwise <- compare_means(
        metric ~ group, data = df,
        method = "wilcox.test"
      ) %>%
        mutate(y.position = max(df$metric)*seq(1.05,1.55,length.out=6))
      
      stat.pairwise$y.position[6] <- stat.pairwise$y.position[4]
      stat.pairwise$y.position[5] <- stat.pairwise$y.position[4]
      stat.pairwise <- stat.pairwise[c(1:3,5),]
      if(input$variablecomp=="Age effect"){
        stat.pairwise <- stat.pairwise[c(1:3,6),]
      }
      
      # Linear mixed effects model
      df_st <- data.frame(values=df$metric,Genotypes=df$group,blocking=as.factor(df$blocking))
      m1 <- mixed(values~Genotypes+(1|blocking),data=df_st)
      # full <- nice(m1)
      # print(full)
      #print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
      (emm1 <- emmeans(m1, "Genotypes"))
      prs <- as.data.frame(pairs(emm1, adjust = "holm"))
      print(prs)
      
      #coefs
      shap_wtne <- shapiro.test(df$metric[df$group==unique(df$group)[1]])
      shap_tgne <- shapiro.test(df$metric[df$group==unique(df$group)[2]])
      shap_wtee <- shapiro.test(df$metric[df$group==unique(df$group)[3]])
      shap_tgee <- shapiro.test(df$metric[df$group==unique(df$group)[4]])
      # annolmm <- as.data.frame(full)
      # annolmm[1,3] <- strsplit(annolmm[1,3],' ')[1]
      # annolmm[1,4] <- as.numeric(annolmm[1,4])
      # annolmm <- rename(annolmm,'contrast'='Effect')
      # annolmm <- rename(annolmm,'F/t.ratio'='F')
      annoprs0 <- prs[c(1:3,5),c(1,4,5,6)]
      if(input$variablecomp3=="Age effect"){
        annoprs0 <- prs[c(1:3,6),c(1,4,5,6)]
      }
      print(annoprs0)
      # annoprs$df <- round(annoprs$df, digits=2)
      # annoprs$t.ratio <- round(annoprs$t.ratio, digits=2)
      # annoprs$p.value <- formatC(annoprs$p.value,digits=2,format = 'g')
      annoprs <- data.frame()[1:7,]
      annoprs$df <- c(formatC(groupannodf,digits=2),round(annoprs0$df, digits=2))
      annoprs$t.ratio <- c(formatC(groupannoF,digits=2),round(annoprs0$t.ratio, digits=2))
      annoprs$p.value <- c(formatC(groupannop,digits=2),formatC(annoprs0$p.value,digits=2,format = 'g'))
      
      annoprs <- rename(as.data.frame(annoprs),'F/t.ratio'='t.ratio')
      # print(annolmm)
      # print(annoprs)
      # anno <- rbind(annolmm,annoprs)
      
      pshap <- c(round(shap_wtne$p.value, digits = 2),
                 round(shap_tgne$p.value, digits = 2),
                 round(shap_wtee$p.value, digits = 2),
                 round(shap_tgee$p.value, digits = 2),
                 '',
                 '',
                 '')
      if(input$variablecomp3=="Age effect"){
        pshap <- c(round(shap_wtne$p.value, digits = 2),
                   round(shap_tgne$p.value, digits = 2),
                   round(shap_wtee$p.value, digits = 2),
                   round(shap_tgee$p.value, digits = 2),
                   '',
                   '',
                   '')
      }
      # metricnm <- rep(input$variablerecmet,5)
      group <- c('WTNE',
                 'TGNE',
                 'WTEE',
                 'TGEE',
                 '',
                 '',
                 ''
      )
      if(input$variablecomp3=="Age effect"){
        group <- c('Wild-Type_2mo',
                   'Wild-Type_6mo',
                   'TgDyrk1A_2mo',
                   'TgDyrk1A_6mo',
                   '',
                   '',
                   ''
        )
      }
      means <- c(
        round(mean(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(mean(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      sds <- c(
        round(sd(df_st$values[df_st$Genotypes==group[1]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[2]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[3]]), digits = 2),
        round(sd(df_st$values[df_st$Genotypes==group[4]]), digits = 2),
        '',
        '',
        ''
      )
      
          contrast <- c(
        'Genotype',
        'Treatment',
        'Genotype:Treatment',
        'WTNE - TGNE',
        'WTNE - WTEE',
        'WTNE - TGEE',
        'TGNE - TGEE'
      )
      if(input$variablecomp3=="Age effect"){
        contrast <- c(
          'Genotype',
          'Age',
          'Genotype:Age',
          'Wild-Type_2mo - Wild-Type_6mo',
          'Wild-Type_2mo - TgDyrk1A_2mo',
          'Wild-Type_2mo - TgDyrk1A_6mo',
          'TgDyrk1A_2mo - TgDyrk1A_6mo'
        )
      }
      # anno <- cbind(metricnm,group,means,sds,pshap,anno)
      anno <- cbind(group,means,sds,pshap,contrast,annoprs)
      
      print(anno)
      # print(xtable(anno, caption = "Statistics for..."), include.rownames = FALSE)
      xtable(anno, caption = "Statistics for...")
    },spacing = 'm')
    
    
    # output$tSNE <- renderPlotly({
    #   tsne <- Rtsne(upData(), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500, check_duplicates=FALSE, num_threads=6)
    #   tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = upDatagroups()[,input$variablegroups])
    #   
    #   print(
    #     ggplotly(
    #       ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) +
    #         #xlim(-.03,.025) +
    #         #ylim(-0.05,0.05) +
    #         # theme_minimal()
    #         theme_classic()
    #       #geom_point(aes_string(color=upDatagroups()[,input$variablegroups],alpha=0.1)) 
    #     )
    #   )
    # })
    # 
    # ######################
    # output$Clustering <- renderPlot({
    #   M <- cor(upData())
    #   corrplot(M) +
    #     theme_classic()
    # })
  })
  
  
  
  
  
  # library(rsconnect)
  # rsconnect::deployApp('./',appName = 'EE_TgDyrk1A')
